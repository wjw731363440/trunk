/*************************************************************************
*  Copyright (C) 2018 by Robert Caulk <rob.caulk@gmail.com>  		 *
*  Copyright (C) 2018 by Bruno Chareyre <bruno.chareyre@hmg.inpg.fr>     *
*									 *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
/* This engine is under active development. Experimental only */

//#define THERMAL
#ifdef THERMAL
#pragma once 

#include<core/PartialEngine.hpp>
#include<core/State.hpp>
#include<pkg/dem/ScGeom.hpp>
#include<pkg/common/Dispatching.hpp>

#ifdef FLOW_ENGINE
//#include<pkg/pfv/FlowEngine.hpp>
#include<lib/triangulation/Tesselation.h>
#include<lib/triangulation/FlowBoundingSphere.hpp>
#include "FlowEngine_FlowEngineT.hpp"
#include<pkg/dem/TesselationWrapper.hpp>
#include<lib/triangulation/Network.hpp>
#endif

// This is how to turn a body thermal without data loss. Should be done in a loop by a single function, ofc.
// Yade [10]: s=sphere((0,0,1),100)
// Yade [11]: s.state.vel=(3,3,3)
// Yade [12]: thState = ThermalState()
// Yade [13]: thState.updateAttrs(s.state.dict())
// Yade [14]: s.state=thState
// Yade [15]: s.state.tmp
//  ->  [15]: 0.0
// Yade [16]: s.state.vel
//  ->  [16]: Vector3(3,3,3)

// Shorter yet strictly equivalent
// Yade [21]: s.state=ThermalState().__setstate__( s.state.__getstate__())


class ThermalState: public State {
	public:
		virtual ~ThermalState();
		// State is declared boost::noncopyable, so copy constructor seems nearly impossible. The solution is to update inherited attributes using python as show in preamble
// 		ThermalState& operator= (const State& source) : State(source) {};//FIXME Thermal.cpp:9:33: error: use of deleted function ‘State& State::operator=(const State&)’

		
		YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(ThermalState,State,"preliminary",
		/*attributes*/
		((Real,temp,0,,"temperature of the body"))
		((bool,oldTempSet,false,,"flag to determine which integration method to use"))
		((Real,tempHold,0,,"holds temperature for 2nd order difference"))
		((Real,oldTemp,0,,"change of temp (for thermal expansion)"))
		((Real,stepFlux,0,,"flux during current step"))
		((Real,capVol,0,,"total overlapping volume"))
		((Real,U,0,,"internal energy of the body"))
		((Real,Cp,0,,"internal energy of the body"))
		((Real,k,0,,"thermal conductivity of the body"))
		((Real,alpha,0,,"coefficient of thermal expansion"))
		((bool,Tcondition,false,,"indicates if particle is assigned dirichlet (constant temp) condition"))
		((int,boundaryId,-1,,"identifies if a particle is associated with constant temperature thrermal boundary condition"))
        	((Real,stabilityCoefficient,0,,"sum of solid and fluid thermal resistivities for use in automatic timestep estimation"))
        	((Real,delRadius,0,,"radius change due to thermal expansion"))
		,
		/* extra initializers */
		,
		/* ctor */ createIndex();
		,
		/* py */
	);
	REGISTER_CLASS_INDEX(ThermalState,State);
};

REGISTER_SERIALIZABLE(ThermalState);

class ThermalEngine : public PartialEngine
{
	public:
		typedef TemplateFlowEngine_FlowEngineT<FlowCellInfo_FlowEngineT,FlowVertexInfo_FlowEngineT> FlowEngineT;
		typedef FlowEngineT::Tesselation					Tesselation;
		typedef FlowEngineT::RTriangulation					RTriangulation;
		typedef FlowEngineT::FiniteCellsIterator				FiniteCellsIterator;
		typedef FlowEngineT::CellHandle						CellHandle;
		typedef FlowEngineT::VertexHandle	VertexHandle;
		typedef std::vector<CellHandle>		VectorCell;
		typedef typename VectorCell::iterator		VCellIterator;

	public:
		Scene* scene;
		bool energySet; //initialize the internal energy of the particles
		FlowEngineT* flow;
        	bool timeStepEstimated;
       		Real thermalDT;
        	bool conductionIter;
        	int elapsedIters;
        	int conductionIterPeriod;
        	Real elapsedTime;
        	bool first;
        	bool runConduction;
        	Real maxTimeStep;
		
		virtual ~ThermalEngine();
		virtual void action();
		void makeThermalState();
		void resetBoundaryFluxSums();
		void setConductionBoundary();
		void thermalExpansion();
		void initializeInternalEnergy();
        	void computeNewPoreTemperatures();
        	void computeNewParticleTemperatures();
		void computeSolidFluidFluxes();
		void computeVertexSphericalArea();
		void computeFlux(CellHandle& cell, const shared_ptr<Body>& b, const double surfaceArea);
		void computeSolidSolidFluxes();
        	void timeStepEstimate();
        	const Real computeCellPressureChangeFromDeltaTemp(CellHandle& cell);
        	const Real computeCellPressureChangeFromDeltaVolume(CellHandle& cell);
        	Real getThermalDT() {return thermalDT;}
        	int getConductionIterPeriod() {return conductionIterPeriod;}
        	Real getMaxTimeStep() {return maxTimeStep;}
		YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(ThermalEngine,PartialEngine,"preliminary",
		/*attributes*/
		((bool,advection,true,,"Activates advection"))
		((bool,conduction,true,,"Activates conduction"))
		((bool,thermoMech,true,,"Activates thermoMech"))
        	((bool,fluidThermoMech,true,,"Activates thermoMech"))
        	((bool,solidThermoMech,true,,"Activates thermoMech"))
		((vector<bool>, bndCondIsTemperature, vector<bool>(6,false),,"defines the type of boundary condition for each side. True if temperature is imposed, False for no heat-flux. Indices can be retrieved with :yref:`FlowEngine::xmin` and friends."))
		((vector<double>, thermalBndCondValue, vector<double>(6,0),,"Imposed value of a boundary condition."))
		((vector<double>, thermalBndFlux, vector<double>(6,0),,"Flux through thermal boundary."))
		((bool,boundarySet,false,,"set false after changing boundary conditions"))
		((bool,useKernMethod,false,,"flag to use Kern method for thermal conductivity area calc"))
        	((Real,fluidBeta,0.0002,,"volumetric temperature coefficient m^3/m^3C, default water, <= 0 deactivates"))
        	((double, fluidK, 0.650,,"Thermal conductivity of the fluid."))
        	((double, tsSafetyFactor, 0.8,,"Allow user to control the timstep estimate with a safety factor. Default 0.8"))
		,
		/* extra initializers */
		,
		/* ctor */
		energySet=false;timeStepEstimated=false;thermalDT=0;elapsedTime=0;first=true;runConduction=false;maxTimeStep=10000;
		,
		/* py */
        	.def("getThermalDT",&ThermalEngine::getThermalDT,"let user check estimated thermalDT .")
        	.def("getConductionIterPeriod",&ThermalEngine::getConductionIterPeriod,"let user check estimated conductionIterPeriod .")
        	.def("getMaxTimeStep",&ThermalEngine::getMaxTimeStep,"let user check estimated maxTimeStep.")
	)
	DECLARE_LOGGER;

	
};
REGISTER_SERIALIZABLE(ThermalEngine);


#endif//THERMAL
	
