package plugins.pikeja.colocalizationstudiomod;

import icy.roi.ROI;
import icy.roi.ROI2D;
import icy.type.point.Point3D;

import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Vector;

import plugins.nchenouard.spot.Spot;


public class apparatedSpots {
	Spot s1;
	Spot s2;
	double distance;

	public apparatedSpots( Spot s1, Spot s2, double distance ) {
		this.s1 = s1;
		this.s2=s2;
		this.distance=distance;
	}

	public static ArrayList<ROI> coloc1 = new ArrayList<ROI>();
	public static ArrayList<ROI> coloc2 = new ArrayList<ROI>();

public static ArrayList<apparatedSpots> appDetectConstruction(Vector<Spot> spots,Vector<Spot> spots2)
{
	ArrayList<apparatedSpots> liste_retour = new ArrayList<apparatedSpots>();
	
	if (spots.isEmpty()){}else{
	int nbdeta=spots.size();
	int nbdetb=spots2.size();
	double x_a,y_a,x_b,y_b,z_a,z_b;
	for (int p2 = 0; p2 < nbdetb; p2++) {
		x_b=spots2.get(p2).mass_center.x;y_b=spots2.get(p2).mass_center.y;z_b=spots2.get(p2).mass_center.z;
		x_a=spots.get(0).mass_center.x;y_a=spots.get(0).mass_center.y;z_a=spots.get(0).mass_center.z;
		double min_dist;
		if ((z_a<0) || (z_b<0)){
			min_dist=Math.sqrt(Math.pow(x_a-x_b, 2)+Math.pow(y_a-y_b, 2));
		}
		else{min_dist=Math.sqrt(Math.pow(x_a-x_b, 2)+Math.pow(y_a-y_b, 2)+Math.pow(z_a-z_b, 2));}
		int indice = 0;
		for (int p = 0; p < nbdeta; p++) {
		x_a=spots.get(p).mass_center.x;y_a=spots.get(p).mass_center.y;z_a=spots.get(p).mass_center.z;
		double temp;
		if ((z_a<0) || (z_b<0)){temp = Math.sqrt(Math.pow(x_a-x_b, 2)+Math.pow(y_a-y_b, 2));}
		else {temp = Math.sqrt(Math.pow(x_a-x_b, 2)+Math.pow(y_a-y_b, 2)+Math.pow(z_a-z_b, 2));}
		if (temp<min_dist)
			{min_dist=temp;indice=p;}
		}
		apparatedSpots aS = new apparatedSpots(spots.get(indice), spots2.get(p2), min_dist);
		liste_retour.add(aS);
		}}
	return liste_retour;
}

//crée la sous liste des spots "colocalisés" par rapport au pourcentage calculé statistiquement	
public static ArrayList<apparatedSpots> appDetectSelect(ArrayList<apparatedSpots> liste_app_detect,int ind_max)
{ArrayList<apparatedSpots> liste_retour = new ArrayList<apparatedSpots>();
//détermination de la distance "max" acceptée
//création d'une liste de toutes les distance
ArrayList<Double> distances = new ArrayList<Double>();
for (apparatedSpots aS:liste_app_detect)
{
	distances.add(aS.distance);
}
Collections.sort(distances);
//détermination de l'indice max qui donne la distance max
//int ind_max = (int)(percentage*distances.size());
if (distances.isEmpty()||ind_max==0){} else{
double distance_max=distances.get(ind_max-1);
//création de la sous liste des apparayedSpots dont la distance est < distance_max
for (apparatedSpots aS:liste_app_detect)
{
	if (aS.distance<=distance_max)
	{
		liste_retour.add(aS);
	}
}	}
return liste_retour;
}

//remplit un tableau de ROIs correspondant aux spots  qui colocalisent à partir de la sélection des spots colocalisés
public static void roiColoc3D( int t, ROI[] detections1,ROI[] detections2,ArrayList<apparatedSpots> spotsColoc)
{			
	for (apparatedSpots aS:spotsColoc)
	{
		Spot s1=aS.s1;
		Spot s2=aS.s2;
		double x1 = s1.mass_center.x;double y1 = s1.mass_center.y;double z1 = s1.mass_center.z;
		double x2 = s2.mass_center.x;double y2 = s2.mass_center.y;double z2 = s2.mass_center.z;
		for (ROI roi1:detections1){			
		Point3D p = massCenters.getMassCenter(roi1);
		double x_1=p.getX();double y_1=p.getY();double z_1=p.getZ();
		boolean time1 = ((roi1.getPosition5D().getT()==t)||(roi1.getBounds5D().isInfiniteT()));
		boolean height1 = (z_1==z1)||(roi1.getBounds5D().isInfiniteZ());
		
		if ((x1==x_1)&(y1==y_1)&time1&height1){coloc1.add(roi1);break;}}
		for (ROI roi2:detections2){			
			Point3D p = massCenters.getMassCenter(roi2);
			double x_2=p.getX();double y_2=p.getY();double z_2 = p.getZ();
			boolean time2 = ((roi2.getPosition5D().getT()==t)||(roi2.getBounds5D().isInfiniteT()));
			boolean height2 = (z_2==z2)||(roi2.getBounds5D().isInfiniteZ());
			if ((x2==x_2)&(y2==y_2)&time2&height2){coloc2.add(roi2);break;}}			
	}				
}
//remplit un tableau de ROIs correspondant aux spots 1 qui colocalisent à partir de la sélection des spots colocalisés
	public static void roiColoc2D( int t, ROI[] detections1,ROI[] detections2,ArrayList<apparatedSpots> spotsColoc)
	{			
		for (apparatedSpots aS:spotsColoc)
		{
			Spot s1=aS.s1;
			Spot s2=aS.s2;
			double x1 = s1.mass_center.x;double y1 = s1.mass_center.y;
			double x2 = s2.mass_center.x;double y2 = s2.mass_center.y;
			for (ROI roi1:detections1){			
			Point2D p = massCenters.getMassCenter2D((ROI2D)roi1);
			double x_1=p.getX();double y_1=p.getY();
			boolean time1 = ((roi1.getPosition5D().getT()==t)||(roi1.getBounds5D().isInfiniteT()));
			if ((x1==x_1)&(y1==y_1)&time1){coloc1.add(roi1);break;}}
			for (ROI roi2:detections2){			
				Point2D p = massCenters.getMassCenter2D((ROI2D) roi2);
				double x_2=p.getX();double y_2=p.getY();
				boolean time2 = ((roi2.getPosition5D().getT()==t)||(roi2.getBounds5D().isInfiniteT()));
				if ((x2==x_2)&(y2==y_2)&time2){coloc2.add(roi2);break;}}			
		}				
	}
	
}