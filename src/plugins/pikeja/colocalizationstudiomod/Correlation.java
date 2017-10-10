package plugins.pikeja.colocalizationstudiomod;

import icy.image.IcyBufferedImage;
import icy.roi.ROI;
import icy.roi.ROI2D;
import icy.sequence.Sequence;
import icy.type.collection.array.Array1DUtil;
import icy.type.rectangle.Rectangle3D;

import java.awt.geom.Rectangle2D;

import plugins.adufour.ezplug.EzVarBoolean;
import plugins.adufour.ezplug.EzVarSequence;

public class Correlation {
	
	public static double[] pearson_TCL(Sequence seq1,
			Sequence seq2, int t,int c1,int c2, ROI roi_t) {
		
		if (roi_t==null){
			double[] p=new double[2];p[0]=0;p[1]=1;
			return p;}
		
		int Z = seq1.getSizeZ();
		
		double[] Pearson = new double[2];
		double mean1 = 0;
		double mean2 = 0;
		double var1 = 0;
		double var2 = 0;
		double m1=0;
		double m2=0;
		double sample = 0d;
		for (int z=0;z<Z;z++){
		IcyBufferedImage img1 = seq1.getImage(t, z);
		IcyBufferedImage img2 = seq2.getImage(t, z);
		Rectangle2D roiBounds1 = roi_t.getBounds5D().toRectangle2D();						
			double[] tab1 = Array1DUtil.arrayToDoubleArray(img1.getDataXY(c1),
					img1.isSignedDataType());
			double[] tab2 = Array1DUtil.arrayToDoubleArray(img2.getDataXY(c2),
					img2.isSignedDataType());

			int sizeX = img1.getSizeX();
			
			int minX = Math.max((int) roiBounds1.getMinX(),(int)img1.getBounds().getMinX())+1;
			int minY =  Math.max((int) roiBounds1.getMinY(),(int)img1.getBounds().getMinY())+1;
			int maxX =  Math.min((int) roiBounds1.getMaxX(),(int)img1.getBounds().getMaxX())-1;
			int maxY =  Math.min((int) roiBounds1.getMaxY(),(int)img1.getBounds().getMaxY())-1;
			

			for (int x = minX; x <= maxX; x++) {
				for (int y = minY; y <= maxY; y++) {
					if (roi_t.contains(x, y,z,t,-1)) {
						int off = (y * sizeX) + x;

						m1 += tab1[off] * tab2[off];
						m2 += Math.pow(tab1[off] * tab2[off],2);
						
						mean1 += tab1[off];
						mean2 += tab2[off];
						
						var1 += Math.pow(tab1[off], 2);
						var2 += Math.pow(tab2[off], 2);
						
						sample++;
					}
				}
			}}

			if (sample>0){
			mean1 /= sample;
			mean2 /= sample;
			m1/=sample;m2/=sample;
			var1 /= sample;
			var2 /= sample;
			
			//calcul des moments du coeff de Pearson sous l'hypothese nulle (mean=0)
			double sigma1=Math.sqrt(var1-Math.pow(mean1,2));double sigma2=Math.sqrt(var2-Math.pow(mean2,2));			
			double m2_R=1/(sample);
			
			//calcul du coefficient de Pearson			
			Pearson[0] = (m1 - mean1*mean2)/(sigma1*sigma2);
			//calcul du coeff de Pearson centré réduit et de son skew
			double R_tilde=Pearson[0]/Math.sqrt(m2_R);			
			Pearson[1]=0.5*(1-ErrorFunction.erf(R_tilde/Math.sqrt(2)));}
		else {Pearson[0]=0;Pearson[1]=0;}		

		return Pearson;
	}


		public static double[] ICCS_compute(Sequence seq1,
				Sequence seq2,int c1,int c2, int t, ROI roi_t) {
			
			if (roi_t==null){
				double[] p=new double[2];
				return p;}
			
			int Z = seq1.getSizeZ();
			double[] ICCS = new double[2];
			double mean1 = 0;
			double mean2 = 0;
			double corr1 = 0,corr2=0,corr_inter=0;
			double sample = 0d;
			
			Rectangle2D roiBounds1 = roi_t.getBounds5D().toRectangle2D();
			IcyBufferedImage img1_ = seq1.getFirstImage();
			int minX = Math.max((int) roiBounds1.getMinX(),(int)img1_.getBounds().getMinX())+1;
			int minY =  Math.max((int) roiBounds1.getMinY(),(int)img1_.getBounds().getMinY())+1;
			int maxX =  Math.min((int) roiBounds1.getMaxX(),(int)img1_.getBounds().getMaxX())-1;
			int maxY =  Math.min((int) roiBounds1.getMaxY(),(int)img1_.getBounds().getMaxY())-1;
			
			for (int z=0;z<Z;z++){
			IcyBufferedImage img1 = seq1.getImage(t, z);
			IcyBufferedImage img2 = seq2.getImage(t, z);
							double[] tab1 = Array1DUtil.arrayToDoubleArray(img1.getDataXY(c1),
						img1.isSignedDataType());
				double[] tab2 = Array1DUtil.arrayToDoubleArray(img2.getDataXY(c2),
						img2.isSignedDataType());

				int sizeX = img1.getSizeX();
				for (int x = minX; x <= maxX; x++) {
					for (int y = minY; y <= maxY; y++) {
						if (roi_t.contains(x, y,z+0.0001,t+0.00001,-1)) {
							int off = (y * sizeX) + x;							
							mean1 += tab1[off];
							mean2 += tab2[off];													
							sample++;
						}
					}
				}}
			
				if (sample>0){
				mean1 /= sample;
				mean2 /= sample;}
				
				
				sample=0;
				//calcul des cross correlations
				for (int z=0;z<Z;z++){
					IcyBufferedImage img1 = seq1.getImage(t, z);
					IcyBufferedImage img2 = seq2.getImage(t, z);
									double[] tab1 = Array1DUtil.arrayToDoubleArray(img1.getDataXY(c1),
								img1.isSignedDataType());
						double[] tab2 = Array1DUtil.arrayToDoubleArray(img2.getDataXY(c2),
								img2.isSignedDataType());
						int sizeX = img1.getSizeX();
				for (int x = minX; x <= maxX; x++) {
					for (int y = minY; y <= maxY; y++) {
						if (roi_t.contains(x, y,z+0.0001,t+0.00001,-1)) {
							int off = (y * sizeX) + x;
							
							corr1 += Math.pow(tab1[off]-mean1,2);
							corr2 += Math.pow(tab2[off]-mean2,2);
							corr_inter += (tab1[off]-mean1)*(tab2[off]-mean2);
						}
					}
				}}
				if (sample>0){
				corr1/=sample;corr2/=sample;corr_inter/=sample;}			
				
				//calcul des fractions P1 et P2 de molecules qui interagissent
				if (corr1>0)
				ICCS[0]=corr_inter/corr1;
				
				if (corr2>0)
					ICCS[1]=corr_inter/corr2;				
											
			return ICCS;
		}

public static double[] MandersCoeff(Sequence seq1,
		Sequence seq2, int t, int c1, int c2, Boolean surface, ROI roi1, ROI roi2,ROI roi_t) {
	
	if (roi_t==null){
		double[] p=new double[4];p[0]=0;p[1]=0;p[2]=1;p[3]=1;
		return p;}
		
	double[] Manders = new double[4];
	Rectangle2D roiBounds1 = roi_t.getBounds5D().toRectangle2D();
	IcyBufferedImage img = seq1.getFirstImage();
	int minX = Math.max((int) roiBounds1.getMinX(),(int)img.getBounds().getMinX())+1;
	int minY =  Math.max((int) roiBounds1.getMinY(),(int)img.getBounds().getMinY())+1;
	int maxX =  Math.min((int) roiBounds1.getMaxX(),(int)img.getBounds().getMaxX())-1;
	int maxY =  Math.min((int) roiBounds1.getMaxY(),(int)img.getBounds().getMaxY())-1;
	
	int Z = seq1.getSizeZ();
	
	double manders1=0;double manders2=0;				
	double manders_surface1=0;double manders_surface2=0;		
	
	double intensite1 = 0;
	double intensite_2_1 = 0;
	
	double intensite2 = 0;
	double intensite_2_2 = 0;
	
	double inter1 = 0;
	double inter2 = 0;
	
	
	double volume_inter=0;
	double volume1 = 0;
	double volume2 = 0;
	double volume_roi = 0;
	
	

	for (int z=0;z<Z;z++){	
	IcyBufferedImage img1 = seq1.getImage(t, z);
	IcyBufferedImage img2 = seq2.getImage(t, z);
	double[] tab1 = Array1DUtil.arrayToDoubleArray(img1.getDataXY(c1),
				img1.isSignedDataType());
	double[] tab2 = Array1DUtil.arrayToDoubleArray(img2.getDataXY(c2),
				img2.isSignedDataType());
	int sizeX = img1.getSizeX();
		
		for (int x = minX; x <= maxX; x++) {			
			for (int y = minY; y <= maxY; y++) {
				if (roi_t.contains(x, y,z,t,-1)){
					volume_roi+=1;
				int off = (y * sizeX) + x;
				if (roi1.contains(x, y,z,t,c1)) {							
					intensite1 += tab1[off];
					intensite_2_1 += Math.pow(tab1[off],2);
					volume1+=1;
					if (roi2.contains(x, y,z,t,c2)) {
						inter1+= tab1[off];
						inter2+= tab2[off];
						volume_inter+=1;								
					}}
				if (roi2.contains(x, y,z,t,c2)){
					intensite2+=tab2[off];
					intensite_2_2 += Math.pow(tab2[off],2);
					volume2+=1;							}
				}
			}
		
	}
	}
		
	if (intensite1 != 0) 
		manders1=inter1/intensite1;
		else
		manders1=0;	
		
		if (intensite2 != 0)
		manders2=inter2/intensite2;
		else
		manders2=0;	
			
		
		if (volume1> 0)
			manders_surface1=volume_inter/volume1;
		else
			manders_surface1=0;
			
		if (volume2 > 0)
			manders_surface2=volume_inter/volume2;
		else
		manders_surface2=0;
		
		if (surface==true){Manders[0]=manders_surface1;Manders[1]=manders_surface2;}
		else{Manders[0]=manders1;Manders[1]=manders2;}
		
		//calcul des p_vallues par TCL sur pixel scrambling
		//int blocksize = block.getValue();
		double p2=volume2/volume_roi;												
		double p1=volume1/volume_roi;
		
		double variance_M1=0,variance_M2=0;
		double variance_M1_surface=0,variance_M2_surface=0;
		double mean_1 = p2;double mean_2 = p1;
		
		if (intensite1>0)
			variance_M1 = (p2*(1-p2)*intensite_2_1)/(Math.pow(intensite1,2));
		
		if (volume1>0)
			variance_M1_surface = p2*(1-p2)/(volume1);
		
		if (intensite2>0)
			variance_M2 = (p1*(1-p1)*intensite_2_2)/(Math.pow(intensite2,2));
		
		if (volume2>0)
			variance_M2_surface = p1*(1-p1)/(volume2);
		
		//sous l'hypothese nulle, M1-mu/sigma est normal-centre-reduit
		double R1=0;
		double R2=0;
		double R1_surface=0;double R2_surface=0;
		
		if (variance_M1>0)
			R1=(manders1-mean_1)/Math.sqrt(variance_M1);
						
		if (variance_M1_surface>0)
			R1_surface=(manders_surface1-mean_1)/Math.sqrt(variance_M1_surface);
						
		if (variance_M2>0)
			R2=(manders2-mean_2)/Math.sqrt(variance_M2);
						
		if (variance_M2_surface>0)
			R2_surface=(manders_surface2-mean_2)/Math.sqrt(variance_M2_surface);
															
		if (surface==true){
			Manders[2]=0.5*(1-ErrorFunction.erf(R1_surface/Math.sqrt(2)));
			if (manders_surface1==1)
				Manders[2]=0;
			
			Manders[3]=0.5*(1-ErrorFunction.erf(R2_surface/Math.sqrt(2)));
			if (manders_surface2==1)
				Manders[3]=0;
		}
		else{	
		Manders[2]=0.5*(1-ErrorFunction.erf(R1/Math.sqrt(2)));				
		if (manders1==1)
			Manders[2]=0;				
		
		Manders[3]=0.5*(1-ErrorFunction.erf(R2/Math.sqrt(2)));				
		if (manders2==1)
			Manders[3]=0;
		}	
		return Manders;
}
}
/*	private double[] ICQ_compute (EzVarSequence sequence1,
EzVarSequence sequence2,int t) {

double[] ICQ = new double[2];

double mean1 = 0;
double mean2 = 0;

Sequence seq1 = sequence1.getValue();
Sequence seq2 = sequence2.getValue();
IcyBufferedImage img1 = seq1.getAllImage().get(t);
IcyBufferedImage img2 = seq2.getAllImage().get(t);

ArrayList<ROI2D> roiArrayList1_t = new ArrayList<ROI2D>();

if (rois.getValue()){
//création des 2 ROIs pour les 2 signaux
//List<ROI2D> roiArrayList1 = seq1.getROI2Ds();
//calcul de la ROI 1 
ArrayList<ROI2D> roiArrayList1 = seq1.getROI2Ds();					
		for (ROI2D roi:roiArrayList1){
			if (roi.getT()==t||roi.getT()==-1){
				roi.setC(0);roi.setZ(0);
				if (roi.getT()==-1){
					ROI2D roiclone=(ROI2D) roi.getCopy();roiclone.setT(t);roiArrayList1_t.add(roiclone);roiclone.delete();}
				else{
					roiArrayList1_t.add(roi);}							
			}	
		}
					
		
		if (roiArrayList1_t.isEmpty())
			roiArrayList1_t.add(new ROI2DRectangle(seq1.getBounds2D()));}
				//roiArrayList1_t.add(new ROI2DRectangle());}						
			
									
else
{roiArrayList1_t.add(new ROI2DRectangle(seq1.getBounds2D()));}

//calcul de l'union des ROIs en forcant le canal à 0
//ROI roiUnion1 = ROIUtil.getUnion(roiArrayList1_t);
ROI2D roiUnion1 = roiArrayList1_t.get(0);
for (ROI2D roi:roiArrayList1_t){						
roiUnion1 = (ROI2D)roiUnion1.getUnion(roi);
roiUnion1.setC(0);roiUnion1.setZ(0);roiUnion1.setT(t);
}
Rectangle2D roiBounds1 = roiUnion1.getBounds5D().toRectangle2D();						

if (roiUnion1 instanceof ROI2D){
ROI2D roi2d = roiUnion1;
double[] tab1 = Array1DUtil.arrayToDoubleArray(img1.getDataXY(0),
		img1.isSignedDataType());
double[] tab2 = Array1DUtil.arrayToDoubleArray(img2.getDataXY(0),
		img2.isSignedDataType());

int sizeX = img1.getSizeX();
int minX = (int) roiBounds1.getMinX();
int minY = (int) roiBounds1.getMinY();
int maxX = (int) roiBounds1.getMaxX();
int maxY = (int) roiBounds1.getMaxY();
double sample = 0d;

//calcul de mean 1 / mean 2
for (int x = minX; x <= maxX; x++) {
	for (int y = minY; y <= maxY; y++) {
		if (roi2d.contains(x, y)) {
			int off = (y * sizeX) + x;							
			mean1 += tab1[off];
			mean2 += tab2[off];
			sample++;
		}
	}
}
if (sample>0){
mean1 /= sample;
mean2 /= sample;}

double icq=0;
//calcul de ICQ
for (int x = minX; x <= maxX; x++) {
	for (int y = minY; y <= maxY; y++) {
		if (roi2d.contains(x, y)) {
			int off = (y * sizeX) + x;							
			if ((tab1[off]-mean1)*(tab2[off]-mean2)>0)
			icq += 1;
		}
	}
}

if (sample>0){
icq/=sample;
icq-=0.5;
ICQ[0]=icq;

//test statistique (sign test) variance sample/4
double sig=0.5/Math.sqrt(sample);
ErrorFunction er = new ErrorFunction();
double tem = ErrorFunction.erf(icq/(sig*Math.sqrt(2)));
ICQ[1]=0.5*(1+tem);}
else {ICQ[0]=0;ICQ[1]=0;}
}

return ICQ;
}*/


