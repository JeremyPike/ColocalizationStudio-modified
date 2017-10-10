package plugins.pikeja.colocalizationstudiomod;

import icy.roi.BooleanMask2D;
import icy.roi.BooleanMask3D;
import icy.roi.ROI;
import icy.roi.ROI3D;
import icy.type.point.Point3D;
import icy.type.point.Point5D;

import java.awt.Point;
import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

import plugins.nchenouard.spot.Spot;

public class Ripley3D {
	
public static double[][]  correlation(int dimZ,ROI roi,List<Spot> spots,List<Spot> spots2,double volume, int nb_a,int nb_b,int N,ArrayList<Double> distance,double ratio_zx) {
		
		double result[][]=new double[N-1][3];
		if (roi==null)
		{return result;}
		
		//results[][1]=K et results[][2]=moyenne distances results[][3]=moyenne distances^2
		double delta_K[] = new double[N-1];	
		double distances_moyennes[]=new double[N-1];
		double distances_2_moyennes[]=new double[N-1];
		double compteur[]=new double[N-1];
		
		//on fait la boucle sur toutes les rois pour incrémenter le calcul de la fonction de Ripley
		Point5D pt = roi.getPosition5D();
		ArrayList<Point3D.Integer> polyg = new ArrayList<Point3D.Integer>();
		//il faut extruder si nécessaire la roi
		if (roi.getDimension()==2){			
			for (int z=0;z<dimZ;z++){				
				BooleanMask2D ma=roi.getBooleanMask2D(z,(int) pt.getT(), (int)pt.getC(), true);
				Point[] tab_point = ma.getContourPoints();		
				for (int i=0;i<tab_point.length;i++){
					Point3D.Integer ptInteger = new Point3D.Integer();
					ptInteger.setX((double)tab_point[i].x);ptInteger.setY((double)tab_point[i].y);
					ptInteger.setZ((double)z);
					polyg.add(ptInteger);
					}
			}
				
		}
		else
		{			
		BooleanMask3D ma = ((ROI3D)roi).getBooleanMask(true);
		Point3D.Integer[] tab_point = ma.getContourPoints();		
		for (int i=0;i<tab_point.length;i++)
			polyg.add(tab_point[i]);
		}
		
		//il faut definir le nb local (par ROI) de detections afin de construire les boucles
		int nbdeta=spots.size();
		int nbdetb=spots2.size();
		double x_a,y_a,x_b,y_b,z_a,z_b;
		for (int p = 0; p < nbdeta; p++) {
			double weight = 1;
			// distance du point a au bord de la ROI (polygon)
			//Point pt1 = new Point(spots.get(p).mass_center.x, spots.get(p).mass_center.y);
			x_a=spots.get(p).mass_center.x;y_a=spots.get(p).mass_center.y;z_a=spots.get(p).mass_center.z;
			
			double d = distance2Polygon(x_a,y_a,z_a, polyg,ratio_zx);
			
			for (int p2 = 0; p2 < nbdetb; p2++) {
				x_b=spots2.get(p2).mass_center.x;y_b=spots2.get(p2).mass_center.y;z_b=spots2.get(p2).mass_center.z;
				double temp = Math.sqrt(Math.pow(x_a-x_b, 2)+Math.pow(y_a-y_b, 2)+Math.pow((z_a-z_b)*ratio_zx, 2));				 
					// Calcul du poids (Ripley)
					weight = 1;
					if( temp>d){
						weight = 0.5+(double)(d/(2*temp));
					}
					        //calcul direct de la fonction de correlation de paire delta K
					for (int l = 1; l < N; l++) {
						if ((temp < distance.get(l))&(temp > distance.get(0))) {											
						delta_K[l-1]+=(1 / weight) * volume
								/ (nb_a * nb_b);
						compteur[l-1]+=1;
						distances_moyennes[l-1]+=temp;
						distances_2_moyennes[l-1]+=Math.pow(temp,2);
						break;
						}
						}					
			}
		}
		
		
		for (int l=0; l<N-1; l++){			
			result[l][0]=delta_K[l];
			if (compteur[l]>0){
			result[l][1]=distances_moyennes[l]/compteur[l];
			result[l][2]=distances_2_moyennes[l]/compteur[l];}
			else
			{result[l][1]=0;result[l][2]=0;}
			}
		return result;		
		}
public static double[] variance_theo(int dimZ,ROI roi,List<Spot> spots, List<Spot> spots2,double volume, int nb_a,int nb_b,int N,ArrayList<Double> distance, int N_h,double ratio_zx){

		
		//la premiere colonne de result[N][1] contient la variance de K, la deuxieme colonne contient quand a elle la variance de delta_K
		
		double result[] = new double[N];
		if (roi==null)
		{return result;}
		
		double[][] D_a = new double[nb_a][nb_a];
		double[][][] A_a = new double[nb_a][nb_a][N];
				
		double[][] h_a = new double[nb_a][N];
		double[] sum_h_a = new double[N];	
		
		//on fait la boucle sur toutes les rois pour incrémenter le calcul de la fonction de Ripley
		//on fait la boucle sur toutes les rois pour incrémenter le calcul de la fonction de Ripley
				Point5D pt = roi.getPosition5D();
				ArrayList<Point3D.Integer> polyg = new ArrayList<Point3D.Integer>();
				//il faut extruder si nécessaire la roi
				if (roi.getDimension()==2){			
					for (int z=0;z<dimZ;z++){				
						BooleanMask2D ma=roi.getBooleanMask2D(z,(int) pt.getT(), (int)pt.getC(), true);
						Point[] tab_point = ma.getContourPoints();		
						for (int i=0;i<tab_point.length;i++){
							Point3D.Integer ptInteger = new Point3D.Integer();
							ptInteger.setX((double)tab_point[i].x);ptInteger.setY((double)tab_point[i].y);
							ptInteger.setZ((double)z);
							polyg.add(ptInteger);
							}
					}
						
				}
				else
				{			
				BooleanMask3D ma = ((ROI3D)roi).getBooleanMask(true);
				Point3D.Integer[] tab_point = ma.getContourPoints();		
				for (int i=0;i<tab_point.length;i++)
					polyg.add(tab_point[i]);
				}
								
				//il faut definir le nb local (par ROI) de detections afin de construire les boucles
				int nbdeta=spots.size();
							

		for (int k = 1; k < N; k++) {
						
			double x_a,y_a,z_a;

			for (int i = 0; i < nbdeta; i++) {
				x_a=spots.get(i).mass_center.x;y_a=spots.get(i).mass_center.y;z_a =spots.get(i).mass_center.z; 				
				double dist = distance2Polygon(x_a,y_a,z_a,polyg,ratio_zx);
				
				
				if ((dist < distance.get(k))&(dist > distance.get(0))) {
					double alpha_=(double)(dist/distance.get(k));
					if (alpha_>0){
	                h_a[i][k] =(-4)*Math.pow(alpha_,3)+6*Math.pow(alpha_,2)-3*alpha_+2+6*Math.pow(alpha_,3)*(alpha_*Math.log(2*alpha_)-Math.log(1+alpha_));}
					else{h_a[i][k]=0;}
				} else
				{h_a[i][k] = 1;}
				
				sum_h_a[k] += h_a[i][k];				    
				double x_b,y_b,z_b;
			
				for (int j = i+1; j < nbdeta; j++) {
					x_b=spots.get(j).mass_center.x;y_b=spots.get(j).mass_center.y;z_b = spots.get(j).mass_center.z;
					double temp = Math.sqrt(Math.pow(x_a-x_b, 2)+Math.pow(y_a-y_b, 2)+Math.pow((z_a-z_b)*ratio_zx, 2));				 
					D_a[i][j] = temp;
					D_a[j][i] = D_a[i][j];
					if (D_a[i][j] < 2 * distance.get(k)) {
						A_a[i][j][k] =(Math.PI/12)*(4*distance.get(k)+D_a[i][j])*Math.pow(2*distance.get(k)-D_a[i][j],2); 
						A_a[j][i][k] = A_a[i][j][k];
					}
										
			}												
		}
		}
		
						
		
		for (int k = 0; k < N; k++) {
			double d3 = Math.pow(distance.get(k),3);
			double e = (4/3)*Math.PI * d3;double e2 = e * e;
		    						
			double temp_A = 0;
				               	        
			for (int p = 0; p < nb_a; p++) {
				for (int n = p; n < nb_a; n++) {
					temp_A = temp_A + A_a[p][n][k] * 2;
				}
			}
			double I3 = (temp_A-(e2/volume)*(nb_a*(nb_a-1)))*nb_b/volume;
			double I4=(e*sum_h_a[k]-nb_a*e2/volume)*nb_b/volume;       
			result[k] = Math.pow(volume/(nb_b * nb_a), 2) * (I3 + I4);			

		}									    		
		return result;

	}
public static double[] variance_theo_delta(int dimZ,ROI roi,List<Spot> spots, List<Spot> spots2,double volume, int nb_a,int nb_b,int N,ArrayList<Double> distance,int N_h,double ratio_zx){

	
	//la premiere colonne de result[N][1] contient la variance de K, la deuxieme colonne contient quand a elle la variance de delta_K
	double result[] = new double[N-1];
	
	double[][] D_a = new double[nb_a][nb_a];
	double[][][] A_a = new double[nb_a][nb_a][N];
	double[][][] A_aa = new double[nb_a][nb_a][N];
	
	double[][] h_a = new double[nb_a][N];				
	//pour balayer tous les points on utilise une variable globale M;
	int M=0;
	
	//on fait la boucle sur toutes les rois pour incrémenter le calcul de la fonction de Ripley
			Point5D pt = roi.getPosition5D();
			ArrayList<Point3D.Integer> polyg = new ArrayList<Point3D.Integer>();
			//il faut extruder si nécessaire la roi
			if (roi.getDimension()==2){			
				for (int z=0;z<dimZ;z++){				
					BooleanMask2D ma=roi.getBooleanMask2D(z,(int) pt.getT(), (int)pt.getC(), true);
					Point[] tab_point = ma.getContourPoints();		
					for (int i=0;i<tab_point.length;i++){
						Point3D.Integer ptInteger = new Point3D.Integer();
						ptInteger.setX((double)tab_point[i].x);ptInteger.setY((double)tab_point[i].y);
						ptInteger.setZ((double)z);
						polyg.add(ptInteger);
						}
				}
					
			}
			else
			{			
			BooleanMask3D ma = ((ROI3D)roi).getBooleanMask(true);
			Point3D.Integer[] tab_point = ma.getContourPoints();		
			for (int i=0;i<tab_point.length;i++)
				polyg.add(tab_point[i]);
			}
			
		//il faut definir le nb local (par ROI) de detections afin de construire les boucles
	int nbdeta=spots.size();										
	for (int k = 1; k < N; k++) {
		
		double d2 = Math.pow(distance.get(k), 2);
		double x_i,y_i,z_i,x_j,y_j,z_j,distance_ij;

		for (int i = 0; i < nbdeta; i++) {
			
			x_i=spots.get(i).mass_center.x;y_i=spots.get(i).mass_center.y;z_i=spots.get(i).mass_center.z;
			double dist = distance2Polygon(x_i,y_i,z_i, polyg,ratio_zx);
			
			
			if ((dist < distance.get(k))&(dist > distance.get(0))) {
				double alpha_=(double)(dist/distance.get(k));
				if (alpha_==0){h_a[i+M][k]=0;}else{
                h_a[i+M][k] =(-4)*Math.pow(alpha_,3)+6*Math.pow(alpha_,2)-3*alpha_+2+6*Math.pow(alpha_,3)*(alpha_*Math.log(2*alpha_)-Math.log(1+alpha_));}
			} else
			{h_a[i+M][k] = 1;}
			
			
			
			for (int j = i+1; j < nbdeta; j++) {
				x_j=spots.get(j).mass_center.x;y_j=spots.get(j).mass_center.y;z_j = spots.get(j).mass_center.z;
				double temp = Math.sqrt(Math.pow(x_i-x_j, 2)+Math.pow(y_i-y_j, 2)+Math.pow((z_i-z_j)*ratio_zx, 2));				 
				D_a[i+M][j+M] = temp;
				D_a[j+M][i+M] = D_a[i+M][j+M];
				if (D_a[i+M][j+M] < 2 * distance.get(k)) {
					A_a[i+M][j+M][k] =(Math.PI/12)*(4*distance.get(k)+D_a[i+M][j+M])*Math.pow(2*distance.get(k)-D_a[i+M][j+M],2); 
					A_a[j+M][i+M][k] = A_a[i+M][j+M][k];
				}										
				
				if (k<N-1){
		            int k2 = k+1; 
					double d2bis=Math.pow(distance.get(k2), 2);					
		                if (D_a[i+M][j+M]<distance.get(k)+distance.get(k2))
		                {
		                	double d1=Math.min(distance.get(k), distance.get(k2));
		                	double d1bis=Math.max(distance.get(k), distance.get(k2));			                                      
		                    if (D_a[i+M][j+M]+d1<d1bis){
		                        A_aa[i+M][j+M][k]=(4/3)*Math.PI*d1*d1*d1;}
		                    else
		                    {
		                    	A_aa[i+M][j+M][k]=(Math.PI/(12*D_a[i+M][j+M]))*(Math.pow(D_a[i+M][j+M],2)+2*D_a[i+M][j+M]*(distance.get(k)+distance.get(k2))-3*Math.pow(distance.get(k)-distance.get(k),2))*(distance.get(k)+distance.get(k2)-Math.pow(D_a[i+M][j+M],2));			                         																															
		                    }			                    
		                    A_aa[j+M][i+M][k]=A_aa[i+M][j+M][k];
		                    
		                }			                			            			            			            								               
				}
		                
		                					
		}												
	}		
	}
	M=M+nbdeta;		                		                 		        		        		        		       								       
							
	
	for (int k = 1; k < N; k++) {
		
		double e2=Math.PI*(4/3)*Math.pow(distance.get(k),3);
		double e1=Math.PI*(4/3)*Math.pow(distance.get(k-1),3);
                
		double temp_A1=0;
        double temp_A2=0;
        double temp_A3=0;
       
        double sum_h_a=0;
        double sum_h_a_bis=0;
        
        for (int p = 0; p < nb_a; p++) {
        	sum_h_a += h_a[p][k];
        	sum_h_a_bis += h_a[p][k-1];
			for (int n = p+1; n < nb_a; n++) {
				temp_A1=temp_A1+A_a[p][n][k-1]*2;
                temp_A2=temp_A2+A_a[p][n][k]*2;	                
                temp_A3=temp_A3+A_aa[p][n][k-1]*4;
			}
		}
		        
        
        
               	       
        double I2 = (temp_A1+temp_A2-temp_A3-(Math.pow(e1,2)/volume+Math.pow(e2,2)/volume-2*e1*e2/volume)*(nb_a*(nb_a-1)))*nb_b/volume;        	      		        
        
        double I1=(e2*sum_h_a-e1*sum_h_a_bis-nb_a*Math.pow(e2-e1,2)/volume)*nb_b/volume;
        	        
        result[k-1] = Math.pow(volume/(nb_b * nb_a), 2) * (I1 + I2) ;	        	

	}			
			
	return result;

}

	public static double distance2Polygon(double x_point, double y_point,double z_point, ArrayList<Point3D.Integer> polyg,double ratio_zx) {
		double dist = Integer.MAX_VALUE;// roiPolyg.getPerimeter();

		int nbc = polyg.size();

		for (int i = 0; i < nbc; i++) {
			Point3D pt1 = polyg.get(i);
			double disttmp = Math.sqrt(Math.pow(x_point-pt1.getX(),2)+Math.pow(y_point-pt1.getY(),2)+Math.pow((z_point-pt1.getZ())*ratio_zx,2));
			dist = Math.min(dist, disttmp);
		}

		return (dist);

	}


}



