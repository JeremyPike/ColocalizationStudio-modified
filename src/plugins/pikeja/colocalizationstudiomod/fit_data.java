package plugins.pikeja.colocalizationstudiomod;

import java.util.ArrayList;

import flanagan.analysis.Regression;
import flanagan.analysis.Stat;

public class fit_data {

	public static double[] coeffs;
	public static double[] residuals;
	//on va mettre les parametres d'initialisation en entrée
	    public static void main(ArrayList<Double> distance,double[] yArray,double area,int nba, double[] start_est){

	            	 
	    	int N=distance.size();
	    	double[] xArray = new double[N];
	    	for (int i=0;i<N;i++)
	    	{
	    		xArray[i]=distance.get(i);
	    	}
	    	
	             // estimates of the standard deviations of y
	             //double[] sdArray = {0.5,0.45,0.55,0.44,0.46,0.51,0.56,0.48,0.5,0.45,0.55,0.44,0.46,0.51,0.56,0.48};


	             // Create instances of the class holding the function, y = a + b.exp(-c.x), evaluation method
	             function_fit f1 = new function_fit();
	             f1.area=area;
	             f1.n_a=nba;	             
	             

	             // assign value to constant b in the function
	             //f1.setB(8.0D);

	             // initial estimates of a and c in y = a + b.exp(-c.x)
	             double[] start = new double[3];
	             

	             start[0] = start_est[0];      // initial estimate of alpha
	             start[1] = start_est[1];      // initial estimate of mu
	             start[2] = start_est[2];      // initial estimate of sigma

	             
	             // initial step sizes for a and c in y = a + b.exp(-c.x)
	             double[] step = new double[3];
	             step[0] = 0.01D;      // initial step size for alpha
	             step[1] = 0.01D;     // initial step size for mu
	             step[2] = 0.01D;     // initial step size for sigma

	             // create an instance of Regression
	             //Regression reg = new Regression(xArray, yArray, sdArray);
	             Regression reg = new Regression(xArray, yArray);

	             // call non-linear regression using default tolerance and maximum iterations and plot display option
	             //alpha compris entre 0 et 1
	             reg.addConstraint(0, -1, 0);
	             reg.addConstraint(0, +1, 1);
	             //reg.addConstraint(0, -1, 0);
	             
	             //mu et sigma >0
	             reg.addConstraint(1, -1, 0);	             
	             reg.addConstraint(2, -1, 0);
	             //mu < r_max	             
	             //reg.simplexPlot(f1, start, step);
	             //reg.simplex(f1, start);
	             //reg.simplex(f1, start, step, 0.0000001, 1000);
	             reg.simplex(f1, start);
	             //int nh=reg.getNiter();
	             coeffs=reg.getCoeff();
	             residuals=reg.getResiduals();
	             
	             
	    }
	    public static void main_manual(ArrayList<Double> distance,double[] yArray,double area,int nba){

       	 
	    	int N=distance.size();
	    	double[] xArray = new double[N];
	    	for (int i=0;i<N;i++)
	    	{
	    		xArray[i]=distance.get(i);
	    	}
	    	
	            

	             // Create instances of the class holding the function, y = a + b.exp(-c.x), evaluation method
	             function_fit f1 = new function_fit();
	             f1.area=area;
	             f1.n_a=nba;	             
	             
	             
	             //triple boucle sur alpha/mu et sigma
	             double alpha_0=0;double delta_alpha=1;
	             double mu_0=0.0;double delta_mu=xArray[xArray.length-1];
	             double sigma_0=0.001;double delta_sigma=delta_mu;
	             	             
	             int N_alpha = 100;int N_mu=100;int N_sigma=100;
	             
	             
	             double alpha_f=alpha_0;double mu_f=mu_0;double sigma_f=sigma_0;
	             
	             double[] f_cac = estimation_function(alpha_0, mu_0, sigma_0, xArray, area, nba);
	             double min_f=estimation_moindre_carre(f_cac,yArray);
	             double[] residus = manualResiduals(f_cac, yArray);
	             //estimation par moindre carré par rapport à Kfit

	             double temp_alpha;double temp_mu;double temp_sigma;
	             for (int i=0;i<N_alpha+1;i++){
	            	 temp_alpha=alpha_0+(i*delta_alpha)/N_alpha;
	            	 for (int j=0;j<N_mu+1;j++){
	            		 temp_mu = mu_0+(j*delta_mu)/N_mu;
	            		 for (int k = 0;k<N_sigma;k++){
	            			 temp_sigma = sigma_0+(k*delta_sigma)/N_sigma;
	            			double[] f_temp=estimation_function(temp_alpha, temp_mu, temp_sigma, xArray, area, nba);
	            			double temp= estimation_moindre_carre(f_temp,yArray);
	            			if (temp<min_f)
	            			{
	            				alpha_f=temp_alpha;mu_f=temp_mu;sigma_f=temp_sigma;
	            				min_f=temp;
	            				residus = manualResiduals(f_temp, yArray);
	            			}
	            			 
	            		 }
	            	 }
	             }
	             double[] coeffs_ = new double[3];
	             coeffs_[0]=alpha_f;coeffs_[1]=mu_f;coeffs_[2]=sigma_f;
	             coeffs=coeffs_;
	             residuals=residus;
	             
	    }
	    
	    public static double[] estimation_function(double temp_alpha,double temp_mu,double temp_sigma,double[] x_array,double area,int n_a){
	    	double[] f=new double[x_array.length];
	    	for (int i=0;i<x_array.length;i++)
	    	{f[i]=temp_alpha*area*0.5/(n_a)*(ErrorFunction.erf((x_array[i] - temp_mu)/(temp_sigma*Math.sqrt(2)))-ErrorFunction.erf((-x_array[i] - temp_mu)/(temp_sigma*Math.sqrt(2))));}	    	
	    	return f;
	    }
	    
	    public static double estimation_moindre_carre(double[] y, double[] y_array){
	    	double t=0;
	    	for (int i=0;i<y_array.length;i++)
	    	{t+=Math.pow(y[i]-y_array[i], 2);}	    	
	    	return t;	    	
	    }

	    public static double[] manualResiduals(double[] y, double[] yArray){	    	
		    double[] temp = new double[y.length];
		    for(int i=0; i<y.length; i++)
		    	temp[i]=yArray[i]-y[i];
		    return temp;
		}
}
