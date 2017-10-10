package plugins.pikeja.colocalizationstudiomod;

import icy.roi.BooleanMask2D;
import icy.roi.ROI;
import icy.roi.ROI2D;

import java.awt.Point;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Vector;

import plugins.nchenouard.spot.Spot;

public class Ripley2D {
	// Calcul de la fonction de correlation (Ripley)
	public static double[][] correlation(
			ROI roi, List<Spot> spots,
			List<Spot> spots2, double volume, int nb_a,
			int nb_b, int N, ArrayList<Double> distance) {

		//ici on peut caster en ROI2D			
				BooleanMask2D ma = ((ROI2D)roi).getAsBooleanMask(true);
				Point[] tab_point = ma.getEdgePoints();
				ArrayList<Point> polyg = new ArrayList<Point>();
				for (int i=0;i<tab_point.length;i++)
					polyg.add(tab_point[i]);
				
				double result[][]=new double[N-1][3];
				//results[][1]=K et results[][2]=moyenne distances results[][3]=moyenne distances^2
				double delta_K[] = new double[N-1];	
				double distances_moyennes[]=new double[N-1];
				double distances_2_moyennes[]=new double[N-1];
				double compteur[]=new double[N-1];
				
		// il faut definir le nb local (par ROI) de detections afin de
			// construire les boucles
			int nbdeta = spots.size();
			int nbdetb = spots2.size();
			double x_a, y_a, x_b, y_b;
			for (int p = 0; p < nbdeta; p++) {
				double weight = 1;
				// distance du point a au bord de la ROI (polygon)
				// Point pt1 = new Point(spots.get(p).mass_center.x,
				// spots.get(p).mass_center.y);
				x_a = spots.get(p).mass_center.x;
				y_a = spots.get(p).mass_center.y;

				double d = distance2Polygon(x_a, y_a, polyg);

				for (int p2 = 0; p2 < nbdetb; p2++) {
					x_b = spots2.get(p2).mass_center.x;
					y_b = spots2.get(p2).mass_center.y;
					double temp = Math.sqrt(Math.pow(x_a - x_b, 2)
							+ Math.pow(y_a - y_b, 2));
					// Calcul du poids (Ripley)
					weight = 1;
					if (temp > d) {
						weight = 1 - (Math.acos(d / temp)) / Math.PI;
					}

					/*
					 * for (int l=0; l<N; l++){ if (temp<distance[l]) { K[l] =
					 * K[l] + (1 / weight) * area / (nb_a * nb_b); } }
					 */
					// calcul direct de la fonction de correlation de paire
					// delta K
					for (int l = 1; l < N; l++) {
						if ((temp < distance.get(l))&(temp > distance.get(0))) {
							delta_K[l - 1] += (1 / weight) * volume
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

	public static double[] variance_theo(ROI roi, List<Spot>spots,
		 List<Spot> spots2, double volume, int nb_a,
			int nb_b, int N, ArrayList<Double> distance, int N_h, double[] results) {

		// la premiere colonne de result[N][1] contient la variance de K, la
		// deuxieme colonne contient quand a elle la variance de delta_K
		double result[] = new double[N];

		double[][] D_a = new double[nb_a][nb_a];
		double[][][] A_a = new double[nb_a][nb_a][N];

		double[][] h_a = new double[nb_a][N];
		double sum_h_a = 0;
		//ici on peut caster en ROI2D			
		BooleanMask2D ma = ((ROI2D)roi).getAsBooleanMask(true);
		Point[] tab_point = ma.getEdgePoints();
		ArrayList<Point> polyg = new ArrayList<Point>();
		for (int i=0;i<tab_point.length;i++)
			polyg.add(tab_point[i]);
					
		int nbdeta = spots.size();

					
			for (int k = 1; k < N; k++) {

				double d2 = Math.pow(distance.get(k), 2);
				double x_a, y_a;

				for (int i = 0; i < nbdeta; i++) {
					x_a = spots.get(i).mass_center.x;
					y_a = spots.get(i).mass_center.y;
					double dist = distance2Polygon(x_a, y_a, polyg);

					if ((dist < distance.get(k))&(dist > distance.get(0))) {
						h_a[i][k] = results[(int) Math.ceil(N_h * dist
								/ distance.get(k))];
					} else {
						h_a[i][k] = 1;
					}

					sum_h_a += h_a[i][k];
					double x_b, y_b;

					for (int j = i + 1; j < nbdeta; j++) {
						x_b = spots.get(j).mass_center.x;
						y_b = spots.get(j).mass_center.y;
						double temp = Math.sqrt(Math.pow(x_a - x_b, 2)
								+ Math.pow(y_a - y_b, 2));

						D_a[i][j] = temp;
						D_a[j][i] = D_a[i][j];
						if (D_a[i][j] < 2 * distance.get(k)) {
							A_a[i][j][k] = 2
									* d2
									* Math.acos(D_a[i][j]
											/ (2 * distance.get(k)))
									- 0.5
									* D_a[i][j]
									* Math.sqrt(4 * d2 - D_a[i][j]
											* D_a[i][j]);
							A_a[j][i][k] = A_a[i][j][k];
						}

					}
				}
			}
			

		for (int k = 0; k < N; k++) {
			double d2 = Math.pow(distance.get(k), 2);

			double e = Math.PI * d2;
			double e2 = e * e;

			double temp_A = 0;

			for (int p = 0; p < nb_a; p++) {
				for (int n = p; n < nb_a; n++) {
					temp_A = temp_A + A_a[p][n][k] * 2;
				}
			}

			double I1 = (temp_A - e2 / volume * (nb_a * (nb_a - 1))) * nb_b
					/ volume;
			double I2 = (e * sum_h_a - nb_a * e2 / volume) * nb_b / volume;
			result[k] = Math.pow(volume / (nb_b * nb_a), 2) * (I1 + I2);

		}
		return result;

	}

	public static double[] variance_theo_delta(ROI roi_t, List<Spot> spots,List<Spot> spots2,double area, int nb_a,int nb_b,int N,ArrayList<Double> distance,int N_h, double[] results){

		
		//la premiere colonne de result[N][1] contient la variance de K, la deuxieme colonne contient quand a elle la variance de delta_K
		double result[] = new double[N-1];
		
		double[][] D_a = new double[nb_a][nb_a];
		double[][][] A_a = new double[nb_a][nb_a][N];
		//double[][][][] A_aa = new double[nb_a][nb_a][N][N];
		double[][][] A_aa = new double[nb_a][nb_a][N];
		
		double[][] h_a = new double[nb_a][N];				
		//pour balayer tous les points on utilise une variable globale M;
		int M=0;		
		//on fait la boucle sur toutes les rois pour incrémenter le calcul de la fonction de Ripley
		
			//on peut caster en ROI2D
		BooleanMask2D ma = ((ROI2D)roi_t).getBooleanMask(true);		
		Point[] liste_p = ma.getContourPoints();
		ArrayList<Point> polyg = new ArrayList<Point>();
		for (int i=0;i<liste_p.length;i++)
			polyg.add(liste_p[i]);
		
		//il faut definir le nb local (par ROI) de detections afin de construire les boucles
		int nbdeta=spots.size();										
		for (int k = 1; k < N; k++) {
			double distancek=distance.get(k);
			double d2 = Math.pow(distancek, 2);
			double x_i,y_i,x_j,y_j,distance_ij;

			for (int i = 0; i < nbdeta; i++) {
				Spot spoti = spots.get(i);
				x_i=spoti.mass_center.x;y_i=spoti.mass_center.y;
				double dist = distance2Polygon(x_i,y_i, polyg);
				
				
		        	if ((dist<distance.get(k))&(dist>distance.get(0))){					
							h_a[i+M][k] = results[(int) Math.ceil(N_h * dist / distancek)];
							} else
						{h_a[i+M][k] = 1;}
		    
				for (int j = i+1; j < nbdeta; j++) {
					Spot spotj = spots.get(j);
					x_j=spotj.mass_center.x;y_j=spotj.mass_center.y;
					distance_ij=Math.sqrt(Math.pow(x_i-x_j,2)+Math.pow(y_i-y_j,2));
					
					D_a[i+M][j+M] = distance_ij;
					D_a[j+M][i+M] = D_a[i+M][j+M];
					if (D_a[i+M][j+M] < 2 * distancek) {
						A_a[i+M][j+M][k] = 2 * d2 * Math.acos(D_a[i+M][j+M] / (2 * distancek)) - 0.5 * D_a[i+M][j+M]
								* Math.sqrt(4 * d2 - D_a[i+M][j+M] * D_a[i+M][j+M] );												
						A_a[j+M][i+M][k] = A_a[i+M][j+M][k];
					}
					
					
					if (k<N-1){
			            int k2 = k+1; 
			            double distancek2=distance.get(k2);
						double d2bis=Math.pow(distancek2, 2);					
			                if (D_a[i+M][j+M]<distancek+distancek2)
			                {
			                	double d1=Math.min(distancek, distancek2);
			                	double d1bis=Math.max(distancek, distancek2);			                                      
			                    if (D_a[i+M][j+M]+d1<d1bis){
			                        A_aa[i+M][j+M][k]=Math.PI*d1*d1;}
			                    else
			                    {
			                    	A_aa[i+M][j+M][k]=d2*Math.acos((Math.pow(D_a[i+M][j+M],2)+d2-d2bis)/(2*D_a[i+M][j+M]*distancek))+d2bis*Math.acos((Math.pow(D_a[i+M][j+M],2)+d2bis-d2)/(2*D_a[i+M][j+M]*distancek2));
			                    	A_aa[i+M][j+M][k]=A_aa[i+M][j+M][k]-0.5*Math.sqrt(((-D_a[i+M][j+M]+distancek+distancek2)*(D_a[i+M][j+M]-distancek+distancek2)*(D_a[i+M][j+M]+distancek-distancek2)*(D_a[i+M][j+M]+distancek+distancek2)));
			                    }			                    
			                    A_aa[j+M][i+M][k]=A_aa[i+M][j+M][k];
			                    
			                }			                			            			            			            								               
					}
			                
			                
					///////////////////////////////////////////////////////////////////
			}												
		}		
		}
		M=M+nbdeta;		                		                 		        		        		        		       								       
										
		for (int k = 1; k < N; k++) {
			
			double e2=Math.PI*Math.pow(distance.get(k),2);
			double e1=Math.PI*Math.pow(distance.get(k-1),2);
	                
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
			        
	        
	        
	               	       
	        double I2 = (temp_A1+temp_A2-temp_A3-(Math.pow(e1,2)/area+Math.pow(e2,2)/area-2*e1*e2/area)*(nb_a*(nb_a-1)))*nb_b/area;        	      		        
	        
	        double I1=(e2*sum_h_a-e1*sum_h_a_bis-nb_a*Math.pow(e2-e1,2)/area)*nb_b/area;
	        	        
	        result[k-1] = Math.pow(area/(nb_b * nb_a), 2) * (I1 + I2) ;	        	
	
		}			
				
		return result;

	}

	public static double distance2Polygon(double x_point, double y_point, ArrayList<Point> polyg) {
		double dist = Integer.MAX_VALUE;// roiPolyg.getPerimeter();

		int nbc = polyg.size();

		for (int i = 0; i < nbc; i++) {
			Point pt1 = polyg.get(i);
			double disttmp = Math.sqrt(Math.pow(x_point-pt1.x,2)+Math.pow(y_point-pt1.y,2));
			dist = Math.min(dist, disttmp);
		}

		return (dist);

	}


}
