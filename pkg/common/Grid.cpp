/*************************************************************************
*  Copyright (C) 2012 by François Kneib   francois.kneib@gmail.com       *
*  Copyright (C) 2012 by Bruno Chareyre   bruno.chareyre@hmg.inpg.fr     *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "Grid.hpp"

//!##################	SHAPES   #####################

GridNode::~GridNode(){}
YADE_PLUGIN((GridNode));

GridConnection::~GridConnection(){}
YADE_PLUGIN((GridConnection));

GridNodeGeom6D::~GridNodeGeom6D(){}
YADE_PLUGIN((GridNodeGeom6D));

ScGridCoGeom::~ScGridCoGeom(){}
YADE_PLUGIN((ScGridCoGeom));

GridCoGridCoGeom::~GridCoGridCoGeom(){}
YADE_PLUGIN((GridCoGridCoGeom));

void GridNode::addConnection(shared_ptr<Body> GC){
	ConnList.push_back(GC);
}

Vector3r GridConnection::getSegment(){
	if (!periodic) return node2->state->pos - node1->state->pos;
 	//else
	const Scene* scene=Omega::instance().getScene().get();
	return node2->state->pos + scene->cell->hSize*cellDist.cast<Real>() - node1->state->pos;
}

Real GridConnection::getLength(){
	return getSegment().norm();
}



//!##################	IGeom Functors   #####################

//!			O-O
bool Ig2_GridNode_GridNode_GridNodeGeom6D::go( const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c)
{	
	//GridConnection* GC = static_cast<GridConnection*>(cm.get());
	bool isNew = !c->geom;
	if (Ig2_Sphere_Sphere_ScGeom::go(cm1,cm2,state1,state2,shift2,force,c)){//the 3 DOFS from ScGeom are updated here
 		if (isNew) {//generate a 6DOF interaction from the 3DOF one generated by Ig2_Sphere_Sphere_ScGeom
			shared_ptr<GridNodeGeom6D> sc (new GridNodeGeom6D());
			*(YADE_PTR_CAST<ScGeom>(sc)) = *(YADE_PTR_CAST<ScGeom>(c->geom));
			c->geom=sc;
		}
		if (updateRotations) YADE_PTR_CAST<GridNodeGeom6D>(c->geom)->precomputeRotations(state1,state2,isNew,creep);
		if(YADE_PTR_CAST<GridNodeGeom6D>(c->geom)->connectionBody){	//test this because the connectionBody may not have been yet initialized.
			YADE_PTR_CAST<GridNodeGeom6D>(c->geom)->connectionBody->state->pos=state1.pos;
		}
		return true;
	}
	else return false;
}

bool Ig2_GridNode_GridNode_GridNodeGeom6D::goReverse( const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c)
{
	return go(cm1,cm2,state2,state1,-shift2,force,c);
}
YADE_PLUGIN((Ig2_GridNode_GridNode_GridNodeGeom6D));

//!			\\//
bool Ig2_GridConnection_GridConnection_GridCoGridCoGeom::go( const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c)
{
	/*FIXME : /!\ Note that this geometry doesn't take care of any unwished duplicated contact or shear force following. /!\*/
	GridConnection* conn1 = YADE_CAST<GridConnection*>(cm1.get());
	GridConnection* conn2 = YADE_CAST<GridConnection*>(cm2.get());
	State* stNode11 = conn1->node1->state.get();
	State* stNode12 = conn1->node2->state.get();
	State* stNode21 = conn2->node1->state.get();
	State* stNode22 = conn2->node2->state.get();
	
	if(conn1->node1==conn2->node1 || conn1->node1==conn2->node2 || conn1->node2==conn2->node1 || conn1->node2==conn2->node2){
		//Two connections share at least one node, so they are contiguous => they must not interact.
		return false;
	}
	//There could be a contact between to connections. Check this now.
	bool isNew = !c->geom;
	Real k,m;
	Vector3r A=stNode11->pos, a=stNode12->pos-A; //"A" is an extremity of conn1, "a" is the connection's segment.
	Vector3r B=stNode21->pos, b=stNode22->pos-B; //"B" is an extremity of conn2, "b" is the connection's segment.
	B+=shift2;//periodicity.
	/* NOW STARTS THE OLD VERSION. IT SHOULD BE REMOVED LATER.
	Vector3r N=a.cross(b);	//"N" is orthogonal to "a" and "b". It means that "N" describes the common plan between a and b.
	if(N.norm()>1e-14){	//If "a" and "b" are colinear, "N==0" and this is a special case.
		Real dist=N.dot(B-A)/(N.norm());	//here "dist" is oriented, so it's sign depends on the orientation of "N" against "AB".
		Vector3r pB=B-dist*(N/(N.norm()));	//"pB" is the projection of the point "B" in the plane defined by his normal vector "N".
		//Now we have pB, so we will compute the intersection of two segments into a plane.
		int b0, b1; //2 base vectors used to compute the segment intersection. For more accuracy and to avoid det==0, don't choose the axis where N is max.
		if(std::abs(N[0])<std::abs(N[1]) || std::abs(N[0])<std::abs(N[2])){b0=0 ; b1=std::abs(N[1])<std::abs(N[2])?1:2;}
		else { b0=1;b1=2;}
		Real det=a[b0]*b[b1]-a[b1]*b[b0];
		if (std::abs(det)>1e-14){
			//Now compute k and m, who are the parameters (relative position on the connections) of the intersection on conn1 ("A" and "a") and conn2 ("B" and "b") respectively.
			k = (b[b1]*(pB[b0]-A[b0])+b[b0]*(A[b1]-pB[b1]))/det;
			m = (a[b0]*(-pB[b1]+A[b1])+a[b1]*(pB[b0]-A[b0]))/det;
			//This is a little bit tricky : if we haven't 0<k,m<1, it means that the intersection is not inside both segments,
			//but the contact can occurs anyway between a connection's extremity and a connection's edge or between two connection's extremity.
			//So the three next lines : don't modify k and m if (0<k,m<1), but modify them otherwise to compute later the right normal and penetrationDepth of the contact.
			k = max(min( k,(Real)1.0),(Real)0.0);
			m = max(min( (A+a*k-B).dot(b)/(pow(b.norm(),2.0)) ,(Real)1.0),(Real)0.0);
			k = max(min( (B+b*m-A).dot(a)/(pow(a.norm(),2.0)) ,(Real)1.0),(Real)0.0);
		}
		else {//should never happen
			k=0;m=0;
			cout<<"Error in Ig2_GridConnection_GridConnection_GridCoGridCoGeom : det=="<<det<<endl;
			cout<<"Details : N="<<N<<" b0="<<b0<<" b1="<<b1<<"  a="<<a<<" b="<<b<<endl;
		}
	}
	else{ //this is a special case for perfectly colinear vectors ("a" and "b")
		Real PA=(A-B).dot(b)/(b.norm()*b.norm()); PA=min((Real)1.0,max((Real)0.0,PA));
		Real Pa=(A+a-B).dot(b)/(b.norm()*b.norm()); Pa=min((Real)1.0,max((Real)0.0,Pa));
		Real PB=(B-A).dot(a)/(a.norm()*a.norm()); PB=min((Real)1.0,max((Real)0.0,PB));
		Real Pb=(B+b-A).dot(a)/(a.norm()*a.norm()); Pb=min((Real)1.0,max((Real)0.0,Pb));
		k=(PB+Pb)/2. ; m=(PA+Pa)/2.;
	} OLD VERSION END*/
	
	/* NOW STARTS THE NEW VERSION */
	Real denom=a.dot(a)*b.dot(b)-pow(a.dot(b),2);
	if(denom!=0){
		k = (a.dot(B-A)*b.dot(b)-a.dot(b)*b.dot(B-A))/denom;
// 		m = (a.dot(b)*a.dot(B-A)-b.dot(B-A)*a.dot(a))/denom; //USELESS BECAUSE DETERMINED FROM k
		k = max(min( k,(Real)1.0),(Real)0.0);
		m = max(min( (A+a*k-B).dot(b)/(pow(b.norm(),2.0)) ,(Real)1.0),(Real)0.0);
		k = max(min( (B+b*m-A).dot(a)/(pow(a.norm(),2.0)) ,(Real)1.0),(Real)0.0);
// 		cout<<"k="<<k<<" m="<<m<<"\n"<<"kc="<<kc<<" mc="<<mc<<"\n\n"<<endl;//}
	}
	else{
		Real PA=(A-B).dot(b)/(b.norm()*b.norm()); PA=min((Real)1.0,max((Real)0.0,PA));
		Real Pa=(A+a-B).dot(b)/(b.norm()*b.norm()); Pa=min((Real)1.0,max((Real)0.0,Pa));
		Real PB=(B-A).dot(a)/(a.norm()*a.norm()); PB=min((Real)1.0,max((Real)0.0,PB));
		Real Pb=(B+b-A).dot(a)/(a.norm()*a.norm()); Pb=min((Real)1.0,max((Real)0.0,Pb));
		k=(PB+Pb)/2. ; m=(PA+Pa)/2.;
	}
	/*NEW VERSION END*/
	
	//Compute the geometry if "penetrationDepth" is positive.
	double penetrationDepth = conn1->radius + conn2->radius - (A+k*a - (B+m*b)).norm();
	shared_ptr<GridCoGridCoGeom> scm;
	if(isNew){
		if(penetrationDepth<0)return false;
		scm=shared_ptr<GridCoGridCoGeom>(new GridCoGridCoGeom());
		c->geom=scm;
	}
	else scm=YADE_PTR_CAST<GridCoGridCoGeom>(c->geom);
	//k and m are used to compute almost everything...
	//Fictious states (spheres) are generated at k or m of each connection, they will handle the contact.
	scm->relPos1=k ; scm->relPos2=m;
	scm->fictiousState1.pos=A + k*a ; scm->fictiousState2.pos=B + m*b;
	scm->radius1 = conn1->radius ; scm->radius2 = conn2->radius;
	scm->fictiousState1.vel = (1-k)*stNode11->vel + k*stNode12->vel;
	scm->fictiousState2.vel = (1-m)*stNode21->vel + m*stNode22->vel;
	Vector3r direction = a/(a.norm());
	scm->fictiousState1.angVel = ((1-k)*stNode11->angVel + k*stNode12->angVel).dot(direction)*direction //twist part : interpolated
	+ a.cross(stNode12->vel - stNode11->vel);// non-twist part : defined from nodes velocities
	direction = b/(b.norm());
	scm->fictiousState2.angVel = ((1-m)*stNode21->angVel + m*stNode22->angVel).dot(direction)*direction //twist part : interpolated
	+ b.cross(stNode22->vel - stNode21->vel);// non-twist part : defined from nodes velocities
	Vector3r normal= scm->fictiousState2.pos - scm->fictiousState1.pos;
	normal/=normal.norm();
	scm->contactPoint = scm->fictiousState1.pos + (scm->radius1-0.5*penetrationDepth)*normal;
	scm->penetrationDepth=penetrationDepth;
	scm->precompute(scm->fictiousState1,scm->fictiousState2,scene,c,normal,isNew,shift2,true);
	return true;
}

bool Ig2_GridConnection_GridConnection_GridCoGridCoGeom::goReverse( const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c)
{
	return go(cm1,cm2,state2,state1,-shift2,force,c);
}
YADE_PLUGIN((Ig2_GridConnection_GridConnection_GridCoGridCoGeom));

//!			O/
bool Ig2_Sphere_GridConnection_ScGridCoGeom::go(	const shared_ptr<Shape>& cm1,
						const shared_ptr<Shape>& cm2,
						const State& state1, const State& state2, const Vector3r& shift2, const bool& force,
						const shared_ptr<Interaction>& c)
{	// Useful variables :
	const State*    sphereSt  = YADE_CAST<const State*>(&state1);
	Sphere*         sphere    = YADE_CAST<Sphere*>(cm1.get());
	GridConnection* gridCo    = YADE_CAST<GridConnection*>(cm2.get());
	GridNode*       gridNo1   = YADE_CAST<GridNode*>(gridCo->node1->shape.get());
	GridNode*       gridNo2   = YADE_CAST<GridNode*>(gridCo->node2->shape.get());
	State*          gridNo1St = YADE_CAST<State*>(gridCo->node1->state.get());
	State*          gridNo2St = YADE_CAST<State*>(gridCo->node2->state.get());
	bool isNew = !c->geom;
	shared_ptr<ScGridCoGeom> scm;
	if (!isNew) scm = YADE_PTR_CAST<ScGridCoGeom>(c->geom);
	else {scm = shared_ptr<ScGridCoGeom>(new ScGridCoGeom());}
	Vector3r segt = gridCo->getSegment();
	Real len = gridCo->getLength();
	Vector3r spherePos = sphereSt->pos - shift2;
	Vector3r branch = spherePos - gridNo1St->pos;
	Vector3r branchN = spherePos - gridNo2St->pos;
	for(int i=0;i<3;i++){
		if(std::abs(branch[i])<1e-14) branch[i]=0.0;
		if(std::abs(branchN[i])<1e-14) branchN[i]=0.0;
	}
	Real relPos = branch.dot(segt)/(len*len);
	if(scm->isDuplicate==2 && scm->trueInt!=c->id2)return true;	//the contact will be deleted into the Law, no need to compute here.
	scm->isDuplicate=0;
	scm->trueInt=-1;
	
	/*
	The 4 conditions below are used to avoid double contact between a sphere and two cylinders, and to follow contact properties when the sphere is sliding along different consecutive GridConnections.
	If none of these conditions are satisfied, the classic contact will be done at the bottom of the Ig2. Else the contact may be copied (if sliding), deleted (if just copied and/or duplicated) and the return statement may be used to abort the Ig2.
	
	The first and the second conditions detect if a sphere's projections is outside the connection. So the contact :
	 - have to be created if the projection is outside all neighbours and not already created.
	 - have to be ignored if the projection is inside at least one neighbour.
	 - if the contact is sliding out to another connection (detected via isNew), mark it as duplicated (it will be ignored by the law and imported (copied) by the new contact).
	 
	 The third and the fourth conditions detect if a sphere's projections is inside the connection. So if the contact is new and :
	  - is before the middle of the connection, we search an old contact that may have slided from one of the previous connections. If we find one, we import it here.
	  - is after the middle of the connection, we search an old contact that may have slided from one of the following connections. If we find one, we import it here.
	 */
	if(relPos<=0){	// if the sphere projection is BEFORE the segment ...
		if(gridNo1->ConnList.size()>1){//	if the node is not an extremity of the Grid (only one connection)
			for(int unsigned i=0;i<gridNo1->ConnList.size();i++){	// ... loop on all the Connections of the same Node ...
				GridConnection* GC = (GridConnection*)gridNo1->ConnList[i]->shape.get();
				if(GC==gridCo)continue;//	self comparison.
				Vector3r segtCandidate1 = GC->node1->state->pos - gridNo1St->pos; // (be sure of the direction of segtPrev to compare relPosPrev.)
				Vector3r segtCandidate2 = GC->node2->state->pos - gridNo1St->pos;
				Vector3r segtPrev = segtCandidate1.norm()>segtCandidate2.norm() ? segtCandidate1:segtCandidate2;
				for(int j=0;j<3;j++){
					if(std::abs(segtPrev[j])<1e-14) segtPrev[j]=0.0;
				}
				Real relPosPrev = (branch.dot(segtPrev))/(segtPrev.norm()*segtPrev.norm());
				// ... and check whether the sphere projection is before the neighbours connections too.
				if(relPosPrev<=0){ //if the sphere projection is outside both the current Connection AND this neighbouring connection, then create the interaction if the neighbour did not already do it before.
					const shared_ptr<Interaction> intr = scene->interactions->find(c->id1,gridNo1->ConnList[i]->getId());
					if(intr && intr->isReal()){
						shared_ptr<ScGridCoGeom> intrGeom=YADE_PTR_CAST<ScGridCoGeom>(intr->geom);
						if(!(intrGeom->isDuplicate==1)){ //skip contact.
							if (isNew) return false;
							else {scm->isDuplicate=1;/*cout<<"Declare "<<c->id1<<"-"<<c->id2<<" as duplicated."<<endl;*/}
						}
					}
				}
				else{//the sphere projection is outside the current Connection but inside the previous neighbour. The contact has to be handled by the Prev GridConnection, not here.
					if (isNew)return false;
					else {
// 						cout<<"The contact "<<c->id1<<"-"<<c->id2<<" may be copied and will be deleted now."<<endl ;
						scm->isDuplicate=1;
						scm->trueInt=-1;
						return true;
					}
				}
			}
		}
	}
	
	//Exactly the same but in the case the sphere projection is AFTER the segment.
	else if(relPos>=1){
		if(gridNo2->ConnList.size()>1){
			for(int unsigned i=0;i<gridNo2->ConnList.size();i++){
				GridConnection* GC = (GridConnection*)gridNo2->ConnList[i]->shape.get();
				if(GC==gridCo)continue;//	self comparison.
				Vector3r segtCandidate1 = GC->node1->state->pos - gridNo2St->pos;
				Vector3r segtCandidate2 = GC->node2->state->pos - gridNo2St->pos;
				Vector3r segtNext = segtCandidate1.norm()>segtCandidate2.norm() ? segtCandidate1:segtCandidate2;
				for(int j=0;j<3;j++){
					if(std::abs(segtNext[j])<1e-14) segtNext[j]=0.0;
				}
				Real relPosNext = (branchN.dot(segtNext))/(segtNext.norm()*segtNext.norm());
				if(relPosNext<=0){ //if the sphere projection is outside both the current Connection AND this neighbouring connection, then create the interaction if the neighbour did not already do it before.
					const shared_ptr<Interaction> intr = scene->interactions->find(c->id1,gridNo2->ConnList[i]->getId());
					if(intr && intr->isReal()){
						shared_ptr<ScGridCoGeom> intrGeom=YADE_PTR_CAST<ScGridCoGeom>(intr->geom);
						if(!(intrGeom->isDuplicate==1)){
							if (isNew) return false;
							else {scm->isDuplicate=1;/*cout<<"Declare "<<c->id1<<"-"<<c->id2<<" as duplicated."<<endl;*/}
						}
 					}
				}
				else{//the sphere projection is outside the current Connection but inside the previous neighbour. The contact has to be handled by the Prev GridConnection, not here.
					if (isNew)return false;
					else {
// 						cout<<"The contact "<<c->id1<<"-"<<c->id2<<" may be copied and will be deleted now."<<endl ;
						scm->isDuplicate=1 ;
						scm->trueInt=-1 ;
						return true;
					}
				}
			}
		}
	}
	
	else if (relPos<=0.5){
		if(gridNo1->ConnList.size()>1){//	if the node is not an extremity of the Grid (only one connection)
			for(int unsigned i=0;i<gridNo1->ConnList.size();i++){	// ... loop on all the Connections of the same Node ...
				GridConnection* GC = (GridConnection*)gridNo1->ConnList[i]->shape.get();
				if(GC==gridCo)continue;//	self comparison.
				Vector3r segtCandidate1 = GC->node1->state->pos - gridNo1St->pos; // (be sure of the direction of segtPrev to compare relPosPrev.)
				Vector3r segtCandidate2 = GC->node2->state->pos - gridNo1St->pos;
				Vector3r segtPrev = segtCandidate1.norm()>segtCandidate2.norm() ? segtCandidate1:segtCandidate2;
				for(int j=0;j<3;j++){
					if(std::abs(segtPrev[j])<1e-14) segtPrev[j]=0.0;
				}
				Real relPosPrev = (branch.dot(segtPrev))/(segtPrev.norm()*segtPrev.norm());
				if(relPosPrev<=0){ //the sphere projection is inside the current Connection and outide this neighbour connection.
					const shared_ptr<Interaction> intr = scene->interactions->find(c->id1,gridNo1->ConnList[i]->getId());
					if( intr && intr->isReal() ){// if an ineraction exist between the sphere and the previous connection, import parameters.
						scm=YADE_PTR_CAST<ScGridCoGeom>(intr->geom);
						if(isNew){
// 							cout<<"Copying contact geom and phys from "<<intr->id1<<"-"<<intr->id2<<" to here ("<<c->id1<<"-"<<c->id2<<")"<<endl;
							c->geom=scm;
							c->phys=intr->phys;
							c->iterMadeReal=intr->iterMadeReal;
						}
						scm->trueInt=c->id2;
						scm->isDuplicate=2;	//command the old contact deletion.
						isNew=0;
						break;
					}
				}
			}
		}
	}
	
	else if (relPos>0.5){
		if(gridNo2->ConnList.size()>1){
			for(int unsigned i=0;i<gridNo2->ConnList.size();i++){
				GridConnection* GC = (GridConnection*)gridNo2->ConnList[i]->shape.get();
				if(GC==gridCo)continue;//	self comparison.
				Vector3r segtCandidate1 = GC->node1->state->pos - gridNo2St->pos;
				Vector3r segtCandidate2 = GC->node2->state->pos - gridNo2St->pos;
				Vector3r segtNext = segtCandidate1.norm()>segtCandidate2.norm() ? segtCandidate1:segtCandidate2;
				for(int j=0;j<3;j++){
					if(std::abs(segtNext[j])<1e-14) segtNext[j]=0.0;
				}
				Real relPosNext = (branchN.dot(segtNext))/(segtNext.norm()*segtNext.norm());
				if(relPosNext<=0){ //the sphere projection is inside the current Connection and outide this neighbour connection.
					const shared_ptr<Interaction> intr = scene->interactions->find(c->id1,gridNo2->ConnList[i]->getId());
					if( intr && intr->isReal() ){// if an ineraction exist between the sphere and the previous connection, import parameters.
						scm=YADE_PTR_CAST<ScGridCoGeom>(intr->geom);
						if(isNew){
// 							cout<<"Copying contact geom and phys from "<<intr->id1<<"-"<<intr->id2<<" to here ("<<c->id1<<"-"<<c->id2<<")"<<endl;
							c->geom=scm;
							c->phys=intr->phys;
							c->iterMadeReal=intr->iterMadeReal;
						}
						scm->trueInt=c->id2;
						scm->isDuplicate=2;	//command the old contact deletion.
						isNew=0;
						break;
					}
				}
			}
		}
	}
	
	relPos=relPos<0?0:relPos;	//min value of relPos : 0
	relPos=relPos>1?1:relPos;	//max value of relPos : 1
	Vector3r fictiousPos=gridNo1St->pos+relPos*segt;
	Vector3r branchF = fictiousPos - spherePos;
 	Real dist = branchF.norm();
	if(isNew && (dist > (sphere->radius + gridCo->radius))) return false;
	
	//	Create the geometry :
	if(isNew) c->geom=scm;
	scm->radius1=sphere->radius;
	scm->radius2=gridCo->radius;
	scm->id3=gridCo->node1->getId();
	scm->id4=gridCo->node2->getId();
	scm->relPos=relPos;
	Vector3r normal=branchF/dist;
	scm->penetrationDepth = sphere->radius+gridCo->radius-dist;
	scm->fictiousState.pos = fictiousPos;
	scm->contactPoint = spherePos + normal*(scm->radius1 - 0.5*scm->penetrationDepth);
	scm->fictiousState.vel = (1-relPos)*gridNo1St->vel + relPos*gridNo2St->vel;
	scm->fictiousState.angVel =
		((1-relPos)*gridNo1St->angVel + relPos*gridNo2St->angVel).dot(segt/len)*segt/len //twist part : interpolated
		+ segt.cross(gridNo2St->vel - gridNo1St->vel);// non-twist part : defined from nodes velocities
	scm->precompute(state1,scm->fictiousState,scene,c,normal,isNew,shift2,true);//use sphere-sphere precompute (with a virtual sphere)
	return true;
}

bool Ig2_Sphere_GridConnection_ScGridCoGeom::goReverse(	const shared_ptr<Shape>& cm1,
								const shared_ptr<Shape>& cm2,
								const State& state1,
								const State& state2,
								const Vector3r& shift2,
								const bool& force,
								const shared_ptr<Interaction>& c)
{
	c->swapOrder();
	return go(cm2,cm1,state2,state1,-shift2,force,c);
}
YADE_PLUGIN((Ig2_Sphere_GridConnection_ScGridCoGeom));


//!##################	Laws   #####################

//!			O/
bool Law2_ScGridCoGeom_FrictPhys_CundallStrack::go(shared_ptr<IGeom>& ig, shared_ptr<IPhys>& ip, Interaction* contact){
	int id1 = contact->getId1(), id2 = contact->getId2();
	ScGridCoGeom* geom= static_cast<ScGridCoGeom*>(ig.get());
	FrictPhys* phys = static_cast<FrictPhys*>(ip.get());
	if(geom->penetrationDepth <0){
		if (neverErase) {
			phys->shearForce = Vector3r::Zero();
			phys->normalForce = Vector3r::Zero();}
		else return false;}
	if (geom->isDuplicate) {
		if (id2!=geom->trueInt) {
			//cerr<<"skip duplicate "<<id1<<" "<<id2<<endl;
			if (geom->isDuplicate==2) return false;
			return true;
		}
	}
	Real& un=geom->penetrationDepth;
	phys->normalForce=phys->kn*std::max(un,(Real) 0)*geom->normal;

	Vector3r& shearForce = geom->rotate(phys->shearForce);
	const Vector3r& shearDisp = geom->shearIncrement();
	shearForce -= phys->ks*shearDisp;
	Real maxFs = phys->normalForce.squaredNorm()*std::pow(phys->tangensOfFrictionAngle,2);

	if (!scene->trackEnergy){//Update force but don't compute energy terms (see below))
		// PFC3d SlipModel, is using friction angle. CoulombCriterion
		if( shearForce.squaredNorm() > maxFs ){
			Real ratio = sqrt(maxFs) / shearForce.norm();
			shearForce *= ratio;}
	} else {
		//almost the same with additional Vector3r instanciated for energy tracing, duplicated block to make sure there is no cost for the instanciation of the vector when traceEnergy==false
		if(shearForce.squaredNorm() > maxFs){
			Real ratio = sqrt(maxFs) / shearForce.norm();
			Vector3r trialForce=shearForce;//store prev force for definition of plastic slip
			//define the plastic work input and increment the total plastic energy dissipated
			shearForce *= ratio;
			Real dissip=((1/phys->ks)*(trialForce-shearForce))/*plastic disp*/ .dot(shearForce)/*active force*/;
			if(dissip>0) scene->energy->add(dissip,"plastDissip",plastDissipIx,/*reset*/false);
		}
		// compute elastic energy as well
		scene->energy->add(0.5*(phys->normalForce.squaredNorm()/phys->kn+phys->shearForce.squaredNorm()/phys->ks),"elastPotential",elastPotentialIx,/*reset at every timestep*/true);
	}
	Vector3r force = -phys->normalForce-shearForce;
	scene->forces.addForce(id1,force);
	scene->forces.addTorque(id1,(geom->radius1-0.5*geom->penetrationDepth)* geom->normal.cross(force));
	Vector3r twist = (geom->radius2-0.5*geom->penetrationDepth)* geom->normal.cross(force);
	scene->forces.addForce(geom->id3,(geom->relPos-1)*force);
	scene->forces.addTorque(geom->id3,(1-geom->relPos)*twist);
	scene->forces.addForce(geom->id4,(-geom->relPos)*force);
	scene->forces.addTorque(geom->id4,geom->relPos*twist);
	return true;
}
YADE_PLUGIN((Law2_ScGridCoGeom_FrictPhys_CundallStrack));


bool Law2_ScGridCoGeom_CohFrictPhys_CundallStrack::go(shared_ptr<IGeom>& ig, shared_ptr<IPhys>& ip, Interaction* contact){
	const int &id1 = contact->getId1();
	const int &id2 = contact->getId2();
	ScGridCoGeom* geom  = YADE_CAST<ScGridCoGeom*> (ig.get());
	CohFrictPhys* phys = YADE_CAST<CohFrictPhys*> (ip.get());
	
	if (geom->isDuplicate) {
		if (id2!=geom->trueInt) {
			//cerr<<"skip duplicate "<<id1<<" "<<id2<<endl;
			if (geom->isDuplicate==2) return false;
			return true;
		}
	}
	
	Vector3r& shearForce    = phys->shearForce;
	if (contact->isFresh(scene) && geom->isDuplicate!=2) shearForce = Vector3r::Zero();
	Real un     = geom->penetrationDepth;
	Real Fn    = phys->kn*(un-phys->unp);
	
	if (phys->fragile && (-Fn)> phys->normalAdhesion) {
		// BREAK due to tension
		return false;
	} else {
		if ((-Fn)> phys->normalAdhesion) {//normal plasticity
			Fn=-phys->normalAdhesion;
			phys->unp = un+phys->normalAdhesion/phys->kn;
			if (phys->unpMax && phys->unp<phys->unpMax)
				return false;
		}
		phys->normalForce = Fn*geom->normal;
		Vector3r& shearForce = geom->rotate(phys->shearForce);
		const Vector3r& dus = geom->shearIncrement();

		//Linear elasticity giving "trial" shear force
		shearForce -= phys->ks*dus;

		Real Fs = phys->shearForce.norm();
		Real maxFs = phys->shearAdhesion;
		if (!phys->cohesionDisablesFriction || maxFs==0)
			maxFs += Fn*phys->tangensOfFrictionAngle;
		maxFs = std::max((Real) 0, maxFs);
		if (Fs  > maxFs) {//Plasticity condition on shear force
			if (phys->fragile && !phys->cohesionBroken) {
				phys->SetBreakingState();
				maxFs = max((Real) 0, Fn*phys->tangensOfFrictionAngle);
			}
			maxFs = maxFs / Fs;
			Vector3r trialForce=shearForce;
			shearForce *= maxFs;
			if (scene->trackEnergy){
				Real dissip=((1/phys->ks)*(trialForce-shearForce))/*plastic disp*/ .dot(shearForce)/*active force*/;
				if(dissip>0) scene->energy->add(dissip,"plastDissip",plastDissipIx,/*reset*/false);}
			if (Fn<0)  phys->normalForce = Vector3r::Zero();//Vector3r::Zero()
		}
		Vector3r force = -phys->normalForce-shearForce;
		scene->forces.addForce(id1,force);
		scene->forces.addTorque(id1,(geom->radius1-0.5*geom->penetrationDepth)* geom->normal.cross(force));
		Vector3r twist = (geom->radius2-0.5*geom->penetrationDepth)* geom->normal.cross(force);
		scene->forces.addForce(geom->id3,(geom->relPos-1)*force);
		scene->forces.addTorque(geom->id3,(1-geom->relPos)*twist);
		scene->forces.addForce(geom->id4,(-geom->relPos)*force);
		scene->forces.addTorque(geom->id4,geom->relPos*twist);
		return true;
	}
}
YADE_PLUGIN((Law2_ScGridCoGeom_CohFrictPhys_CundallStrack));

bool Law2_GridCoGridCoGeom_FrictPhys_CundallStrack::go(shared_ptr<IGeom>& ig, shared_ptr<IPhys>& ip, Interaction* contact){
	int id1 = contact->getId1(), id2 = contact->getId2();
	id_t id11 = (static_cast<GridConnection*>((&Body::byId(id1)->shape)->get()))->node1->getId();
	id_t id12 = (static_cast<GridConnection*>((&Body::byId(id1)->shape)->get()))->node2->getId();
	id_t id21 = (static_cast<GridConnection*>((&Body::byId(id2)->shape)->get()))->node1->getId();
	id_t id22 = (static_cast<GridConnection*>((&Body::byId(id2)->shape)->get()))->node2->getId();
	GridCoGridCoGeom*    geom= static_cast<GridCoGridCoGeom*>(ig.get());
	FrictPhys* phys = static_cast<FrictPhys*>(ip.get());
	if(geom->penetrationDepth <0){
		if (neverErase) {
			phys->shearForce = Vector3r::Zero();
			phys->normalForce = Vector3r::Zero();}
		else return false;}
	Real& un=geom->penetrationDepth;
	phys->normalForce=phys->kn*std::max(un,(Real) 0)*geom->normal;

	Vector3r& shearForce = geom->rotate(phys->shearForce);
	const Vector3r& shearDisp = geom->shearIncrement();
	shearForce -= phys->ks*shearDisp;
	Real maxFs = phys->normalForce.squaredNorm()*std::pow(phys->tangensOfFrictionAngle,2);

	if (!scene->trackEnergy  && !traceEnergy){//Update force but don't compute energy terms (see below))
		// PFC3d SlipModel, is using friction angle. CoulombCriterion
		if( shearForce.squaredNorm() > maxFs ){
			Real ratio = sqrt(maxFs) / shearForce.norm();
			shearForce *= ratio;}
	} else {
		//almost the same with additional Vector3r instatinated for energy tracing, 
		//duplicated block to make sure there is no cost for the instanciation of the vector when traceEnergy==false
		if(shearForce.squaredNorm() > maxFs){
			Real ratio = sqrt(maxFs) / shearForce.norm();
			Vector3r trialForce=shearForce;//store prev force for definition of plastic slip
			//define the plastic work input and increment the total plastic energy dissipated
			shearForce *= ratio;
			Real dissip=((1/phys->ks)*(trialForce-shearForce))/*plastic disp*/ .dot(shearForce)/*active force*/;
			if (traceEnergy) plasticDissipation += dissip;
			else if(dissip>0) scene->energy->add(dissip,"plastDissip",plastDissipIx,/*reset*/false);
		}
		// compute elastic energy as well
		scene->energy->add(0.5*(phys->normalForce.squaredNorm()/phys->kn+phys->shearForce.squaredNorm()/phys->ks),"elastPotential",elastPotentialIx,/*reset at every timestep*/true);
	}
	Vector3r force = -phys->normalForce-shearForce;
	Vector3r torque1 = (geom->radius1-0.5*geom->penetrationDepth)* geom->normal.cross(force);
	Vector3r torque2 = (geom->radius2-0.5*geom->penetrationDepth)* geom->normal.cross(force);
	
	scene->forces.addForce(id11,(1-geom->relPos1)*force);
	scene->forces.addForce(id12,geom->relPos1*force);
	scene->forces.addForce(id21,-(1-geom->relPos2)*force);
	scene->forces.addForce(id22,-geom->relPos2*force);
	
	scene->forces.addTorque(id11,(1-geom->relPos1)*torque1);
	scene->forces.addTorque(id12,geom->relPos1*torque1);
	scene->forces.addTorque(id21,(1-geom->relPos2)*torque2);
	scene->forces.addTorque(id22,geom->relPos2*torque2);
	return true;
}
YADE_PLUGIN((Law2_GridCoGridCoGeom_FrictPhys_CundallStrack));
//!##################	Bounds   #####################

void Bo1_GridConnection_Aabb::go(const shared_ptr<Shape>& cm, shared_ptr<Bound>& bv, const Se3r& se3, const Body* b){
	GridConnection* GC = static_cast<GridConnection*>(cm.get());
	if(!bv){ bv=shared_ptr<Bound>(new Aabb); }
	Aabb* aabb=static_cast<Aabb*>(bv.get());
	Vector3r O = YADE_CAST<State*>(GC->node1->state.get())->pos;
	Vector3r O2 = YADE_CAST<State*>(GC->node2->state.get())->pos;
 	if(!scene->isPeriodic){
		for (int k=0;k<3;k++){
			aabb->min[k]=min(O[k],O2[k])-GC->radius;
			aabb->max[k]=max(O[k],O2[k])+GC->radius;
		}
		return;
 	}
 	else{
		O = scene->cell->unshearPt(O);
		O2 = scene->cell->unshearPt(O2);
		O2 = O2 + scene->cell->hSize*GC->cellDist.cast<Real>();
		for (int k=0;k<3;k++){
			aabb->min[k]=min(O[k],O2[k])-GC->radius;
			aabb->max[k]=max(O[k],O2[k])+GC->radius;
		}
	}
}

YADE_PLUGIN((Bo1_GridConnection_Aabb));
