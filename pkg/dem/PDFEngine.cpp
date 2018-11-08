template<class Phys>
void PDFEngine::getStressSpectrums(vector<vector<vector<Matrix3r> > > &spects, vector<Vector3r Phys::* > const& origins, int nTheta, int nPhi)
{
	const shared_ptr<Scene>& scene=Omega::instance().getScene();
	
	// Resize vectors
	spects.resize(origins.size());
	
	for(uint l(0);l<spects.size();l++) {
		spects[l].resize(nTheta);
		for(int i(0);i<nTheta;i++) {
			spects[l][i].resize(nPhi);			
			for(int j(0);j<nPhi;j++)
				spects[l][i][j] = Matrix3r::Zero();
		}
	}
	
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions) {
		if(!I->isReal()) continue;
		ScGeom* geom=YADE_CAST<ScGeom*>(I->geom.get());
		Phys* phys=YADE_CAST<Phys*>(I->phys.get());
		
		if(phys && geom)
		{
			Real r = geom->refR1 + geom->refR2 - geom->penetrationDepth;
			Vector3r l = r/scene->cell->getVolume()*geom->normal;
			
			// Convention: Oxy: phi=0, Oxz: theta=pi/2, Oyz: phi=pi/2
			Real theta 	= acos(geom->normal.y()); //[0;pi]
			Real phi 	= atan2(geom->normal.z(), geom->normal.x()) + Mathr::PI; //[-pi;pi] => [0; 2pi] for index calculation
			
			Real dTheta = (Mathr::PI / nTheta);
			Real dPhi	= (2.*Mathr::PI / nPhi);
			Real dS = sin(theta) * dTheta * dPhi;
			l = l / dS;
			
			// Calculate indexes
			int idT1 = ((int)round((theta) / dTheta)) % nTheta;
			int idP1 = ((int)round((phi) / dPhi)) % nPhi;
			int idT2 = ((int)round((Mathr::PI - theta) / dTheta)) % nTheta;
			int idP2 = ((int)round((Mathr::PI + phi) / dPhi)) % nPhi;
			
			// Compute stress
			for(uint j(0);j<origins.size();j++) {
				spects[j][idT1][idP1] += phys->*(origins[j])*l.transpose();
				spects[j][idT2][idP2] += phys->*(origins[j])*l.transpose();
			}
		}
	}
}

template<class Phys, class Geom, class datatype>
void PDFEngine::getSimpleSpectrum(vector<vector<vector<datatype> > > & spectrum, vector<PDFEngine::genPtr<Phys, Geom, datatype> > const&, int, int)
{
	
}

void PDFEngine::writeToFile(vector<vector<vector<Matrix3r> > > const& stress_data, vector<string> const& stress_meta,
		vector<vector<vector<Matrix3r> > > const& mat_data, vector<string> const& mat_meta,
		vector<vector<vector<Vector3r> > > const& vector_data, vector<string> const& vec_meta,
		vector<vector<vector<Real> > > const& scalar_data, vector<string> const& scalar_meta)
{
	string const mat_suffix[] = {"xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz"};
	string const vec_suffix[] = {"x","y","z"};
	
	FILE* fid = fopen(filename.c_str(), "a");
	
	if(fid) {
		if(firstRun) {
			fprintf(fid, "# ");
			// Header for stress
			for(uint l(0);l<stress_meta.size();l++)
				for(int i(0);i<numDiscretizeAngleTheta;i++)
					for(int j(0);j<numDiscretizeAnglePhi;j++)
						for(int k(0);k<9;k++)
							fprintf(fid,"%s_%s(%f,%f)\t", stress_meta[l].c_str(), mat_suffix[k].c_str(), i*Mathr::PI/numDiscretizeAngleTheta, j*Mathr::PI*2./numDiscretizeAnglePhi - Mathr::PI);
						
			// Header for matrix
			for(uint l(0);l<mat_meta.size();l++)
				for(int i(0);i<numDiscretizeAngleTheta;i++)
					for(int j(0);j<numDiscretizeAnglePhi;j++)
						for(int k(0);k<9;k++)
							fprintf(fid,"%s_%s(%f,%f)\t", mat_meta[l].c_str(), mat_suffix[k].c_str(), i*Mathr::PI/numDiscretizeAngleTheta, j*Mathr::PI*2./numDiscretizeAnglePhi - Mathr::PI);
						
			// Vector
			for(uint l(0);l<vec_meta.size();l++)
				for(int i(0);i<numDiscretizeAngleTheta;i++)
					for(int j(0);j<numDiscretizeAnglePhi;j++)
						for(int k(0);k<3;k++)
							fprintf(fid,"%s_%s(%f,%f)\t", vec_meta[l].c_str(), vec_suffix[k].c_str(), i*Mathr::PI/numDiscretizeAngleTheta, j*Mathr::PI*2./numDiscretizeAnglePhi - Mathr::PI);
			
			// Scalar
			for(uint l(0);l<scalar_meta.size();l++)
				for(int i(0);i<numDiscretizeAngleTheta;i++)
					for(int j(0);j<numDiscretizeAnglePhi;j++)
						fprintf(fid,"%s(%f,%f)\t", scalar_meta[l].c_str(), i*Mathr::PI/numDiscretizeAngleTheta, j*Mathr::PI*2./numDiscretizeAnglePhi - Mathr::PI);
			fprintf(fid, "\n");
			firstRun = false;
		}

		for(uint l(0);l<stress_data.size();l++)
			for(int i(0);i<numDiscretizeAngleTheta;i++)
				for(int j(0);j<numDiscretizeAnglePhi;j++)
					for(int I(0);I<3;I++)
						for(int J(0);J<3;J++)
							fprintf(fid,"%f\t",stress_data[l][i][j](I,J));
						
		for(uint l(0);l<mat_data.size();l++)
			for(int i(0);i<numDiscretizeAngleTheta;i++)
				for(int j(0);j<numDiscretizeAnglePhi;j++)
					for(int I(0);I<3;I++)
						for(int J(0);J<3;J++)
							fprintf(fid,"%f\t",mat_data[l][i][j](I,J));
						
		for(uint l(0);l<vector_data.size();l++)
			for(int i(0);i<numDiscretizeAngleTheta;i++)
				for(int j(0);j<numDiscretizeAnglePhi;j++)
					for(int I(0);I<3;I++)
						fprintf(fid,"%f\t",vector_data[l][i][j](I));
					
		for(uint l(0);l<scalar_data.size();l++)
			for(int i(0);i<numDiscretizeAngleTheta;i++)
				for(int j(0);j<numDiscretizeAnglePhi;j++)
					fprintf(fid,"%f\t",scalar_data[l][i][j]);

		fprintf(fid, "\n");		
		fclose(fid);
	}
}

CREATE_LOGGER(PDFEngine);
