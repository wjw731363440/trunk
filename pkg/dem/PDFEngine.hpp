
// Should be template class. Ugly but boost::python doesnt like that...
class PDFEngine : public PeriodicEngine {

public:
	
	template<class Phys>
	static void getStressSpectrums(vector<vector<vector<Matrix3r> > > &spects, vector<Vector3r Phys::* > const&, int,int);
	
	template<class Phys, class Geom, class datatype> using genPtr = datatype *(Phys*, Geom*); 
	
	template<class Phys, class Geom, class datatype> 
	static void getSimpleSpectrum(vector<vector<vector<datatype> > > & spectrum, vector<genPtr<Phys, Geom, datatype> > const&, int, int); 
	
	void writeToFile(vector<vector<vector<Matrix3r> > > const& stress_data, vector<string> const&,
		vector<vector<vector<Matrix3r> > > const& mat_data, vector<string> const&,
		vector<vector<vector<Vector3r> > > const& vector_data, vector<string> const&,
		vector<vector<vector<Real> > > const& scalar_data, vector<string> const&);
	virtual void action() {};
	
	YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(PDFEngine, PeriodicEngine,
		"Abstract class for spectrums calculation",
		((int, numDiscretizeAngleTheta, 20,,"Number of sector for theta-angle"))
		((int, numDiscretizeAnglePhi, 40,,"Number of sector for phi-angle"))
		//((Real, discretizeRadius, 0.1,,"d/a interval size"))
		((string, filename, "", , "Filename"))
		((bool, firstRun, true, (Attr::hidden | Attr::readonly), ""))
		,,
		//.def("getSpectrums", &LubricationDPFEngine::PyGetSpectrums,(py::arg("nPhi")=40, py::arg("nTheta")=20), "Get Stress spectrums").staticmethod("getSpectrums")
	);
	DECLARE_LOGGER;
};

REGISTER_SERIALIZABLE(PDFEngine);
