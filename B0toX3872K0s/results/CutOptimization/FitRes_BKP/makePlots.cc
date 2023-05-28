{

	double Xmin = -0.5, Dx = 0.05, X;
	int Nx = 14;
	double mRmin = 0.5, Dm = 0.05, M;
	int Nm = 5;

	for(int kx = 0; kx < Nx; kx++){
		X = Xmin + kx*Dx;
		for (int km = 0; km < Nm; km++){
			M = mRmin + km*Dm;
			std::cout << X << std::endl;
			FitPlot(X, M);	
			CorrMtx(X, M);	

		}	
	}


}
