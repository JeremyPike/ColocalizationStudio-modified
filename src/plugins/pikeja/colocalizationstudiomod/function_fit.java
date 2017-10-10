package plugins.pikeja.colocalizationstudiomod;

import flanagan.analysis.RegressionFunction;
import flanagan.analysis.Stat;

class function_fit implements RegressionFunction{

public double area=87000;
public int n_a=212;  
    @Override
	public double function(double[] p, double[] x){
    	double y=p[0]*area*0.5/(n_a)*(Stat.erf((x[0] - p[1])/(p[2]*Math.sqrt(2)))-Stat.erf((-x[0] - p[1])/(p[2]*Math.sqrt(2))));
         return y;
    }
}

