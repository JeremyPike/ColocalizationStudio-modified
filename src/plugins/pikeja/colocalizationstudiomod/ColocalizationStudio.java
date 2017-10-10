package plugins.pikeja.colocalizationstudiomod;

import java.awt.Point;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.List;
import java.util.Vector;

import javax.swing.JSeparator;
import javax.vecmath.Point3d;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import flanagan.analysis.Stat;
import icy.file.FileUtil;
import icy.gui.frame.progress.AnnounceFrame;
import icy.main.Icy;
import icy.roi.BooleanMask2D;
import icy.roi.BooleanMask3D;
import icy.roi.ROI;
import icy.roi.ROI2D;
import icy.sequence.Sequence;
import icy.sequence.SequenceEvent;
import icy.sequence.SequenceListener;
import icy.system.profile.Chronometer;
import icy.type.point.Point5D;
import icy.util.StringUtil;
import icy.util.XLSUtil;
import jxl.write.WritableSheet;
import jxl.write.WritableWorkbook;
import jxl.write.WriteException;
import plugins.adufour.ezplug.EzButton;
import plugins.adufour.ezplug.EzGroup;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzVarBoolean;
import plugins.adufour.ezplug.EzVarDouble;
import plugins.adufour.ezplug.EzVarEnum;
import plugins.adufour.ezplug.EzVarFile;
import plugins.adufour.ezplug.EzVarInteger;
import plugins.adufour.ezplug.EzVarSequence;
import plugins.adufour.ezplug.EzVarSwimmingObject;
import plugins.adufour.ezplug.EzVarText;
import plugins.adufour.quickhull.QuickHull2D;
import plugins.adufour.quickhull.QuickHull3D;
import plugins.adufour.roi.mesh.polygon.ROI3DPolygonalMesh;
import plugins.kernel.roi.roi2d.ROI2DArea;
import plugins.kernel.roi.roi2d.ROI2DPolygon;
import plugins.kernel.roi.roi2d.ROI2DRectangle;
import plugins.kernel.roi.roi3d.ROI3DArea;
import plugins.nchenouard.spot.DetectionResult;
import plugins.nchenouard.spot.Spot;

// Colocalisation with Ripley function K
// Significant 

public class ColocalizationStudio extends EzPlug implements ActionListener, SequenceListener {

	EzButton hull_b = new EzButton("Define ROI with detections 2 Convex hull ?", this);
	EzButton setTrackSetA = new EzButton("Select ROIs manually ?", this);
	EzButton reset = new EzButton("Reset ROIs ", this);
	Boolean hull = new Boolean(false);
	ROI[] coloc_roi_1 = null;
	ROI[] coloc_roi_2 = null;
	ArrayList<ROI> hull_roi = new ArrayList<ROI>();

	ArrayList<ROI> list_roi = new ArrayList<ROI>();
	boolean no_roi = false;

	private enum ColocMethod {
		CORRELATION, OBJECT
	}

	private EzVarEnum<ColocMethod> method = new EzVarEnum<ColocMethod>("Method", ColocMethod.values(),
			ColocMethod.CORRELATION);

	EzVarSequence sequence1 = new EzVarSequence("Sequence 1");
	EzVarSequence sequence2 = new EzVarSequence("Sequence 2");
	EzVarInteger channel1 = new EzVarInteger("Channel 1");
	EzVarInteger channel2 = new EzVarInteger("Channel 2");

	EzVarSwimmingObject<DetectionResult> detections1 = new EzVarSwimmingObject<DetectionResult>("Detections1");
	EzVarSwimmingObject<DetectionResult> detections2 = new EzVarSwimmingObject<DetectionResult>("Detections2");

	// ROIs (spots) qui colocalisent en sortie
	EzVarBoolean export_colocalized_rois = new EzVarBoolean("Export colocalized detections (ROIs)", false);

	// private EzVarBoolean rois = new EzVarBoolean("Restrict the analysis to
	// ROIs of sequence 1", false) ;

	private EzVarDouble pearson_coeff = new EzVarDouble("Pearson R", 1, 0, 1, 0.1);
	private EzVarText pearson_coeff_label = new EzVarText("Pearson R", "   N/A   ");

	// private EzVarDouble blocksize = new EzVarDouble("Block size for the
	// randomization", 1, 0, 10, 1);
	private EzVarDouble pvalue_pearson = new EzVarDouble("p value Pearson (Randomization)", 1, 0, 1, 0.1);
	private EzVarText pvalue_pearson_label = new EzVarText("p value", "   N/A   ");

	private EzVarBoolean surface = new EzVarBoolean("Surface", false);
	private EzVarDouble M1 = new EzVarDouble("Manders M1", 1, 0, 1, 0.1);
	private EzVarText M1_label = new EzVarText("Manders M1", "   N/A   ");
	private EzVarDouble pvalue_M1 = new EzVarDouble("p value M1 (Randomization)", 1, 0, 1, 0.1);
	private EzVarText pvalue_M1_label = new EzVarText("p value M1", "   N/A   ");

	private EzVarDouble M2 = new EzVarDouble("Manders M2", 1, 0, 1, 0.1);
	private EzVarText M2_label = new EzVarText("Manders M2", "   N/A   ");
	private EzVarDouble pvalue_M2 = new EzVarDouble("p value M2 (Randomization)", 1, 0, 1, 0.1);
	private EzVarText pvalue_M2_label = new EzVarText("p value M2", "   N/A   ");

	/*
	 * private EzVarDouble ICQ=new EzVarDouble("ICQ", 0, -0.5, 0.5, 0.01);
	 * private EzVarText ICQ_label = new EzVarText("ICQ");
	 * 
	 * private EzVarDouble pvalue_ICQ = new EzVarDouble("p value (Sign test)",
	 * 1, 0, 1, 0.1); private EzVarText pvalue_ICQ_label = new
	 * EzVarText("p value");
	 */

	private EzVarDouble ICCS1 = new EzVarDouble("Cross-Correlation 1", 0, 0, 1, 0.01);
	private EzVarText ICCS1_label = new EzVarText("Cross-Correlation 1", "   N/A   ");
	private EzVarDouble ICCS2 = new EzVarDouble("Cross-Correlation 2", 0, 0, 1, 0.01);
	private EzVarText ICCS2_label = new EzVarText("Cross-Correlation 2", "   N/A   ");

	private EzVarDouble max_radius;
	// private EzVarBoolean non_param;
	private EzVarBoolean manual;

	private EzVarDouble number_fit;
	private EzVarText number_fit_label = new EzVarText("Number of colocalized spots 2", "   N/A   ");

	private EzVarDouble alpha_fit;
	private EzVarText alpha_fit_label = new EzVarText("Percentage of colocalized spots 2", "   N/A   ");
	private EzVarDouble mu_fit;
	private EzVarText mu_fit_label = new EzVarText("Mean colocalization distance (pixels)", "   N/A   ");
	private EzVarDouble sigma_fit;
	private EzVarText sigma_fit_label = new EzVarText("Std. dev. of the colocalization distance (pixels)", "   N/A   ");

	private EzVarDouble sup_K;
	private EzVarText sup_K_label = new EzVarText("Max. of the Ripley's K function", "   N/A   ");
	private EzVarDouble p_value_K;
	private EzVarText p_value_K_label = new EzVarText("p-value", "   N/A   ");

	private EzVarBoolean graph = new EzVarBoolean("Plot the the Ripley's K function (fit check)", false);
	// variables globales pour le graphe
	private ArrayList<Double> distance_fit = new ArrayList<Double>();
	private ArrayList<Double> K = new ArrayList<Double>();
	private ArrayList<Double> K_fit = new ArrayList<Double>();
	ChartPanel chartpanel;

	protected EzVarBoolean exportExcel = new EzVarBoolean("Export to Excel", false);
	protected EzVarFile exportExcelFile = new EzVarFile("Excel file", "");

	double[] results = null;
	int N_h;
	// pour le test stat
	double mindist;

	private void betaCorrection(double pas, int nbN) {

		Double valN = 1 / (1 / (double) nbN);
		N_h = valN.intValue() + 1; // / ATTENTION AUX BORNES

		double[] alpha = new double[N_h + 1];
		results = new double[N_h + 1];

		for (int i = 0; i < results.length; i++) {
			results[i] = 0;
			alpha[i] = 0.0d;

		}

		for (int i = 1; i < results.length; i++) {
			alpha[i] = (i / (double) N_h);
		}
		for (int i = 0; i < results.length; i++) {
			double h, j;
			for (h = alpha[i] + (pas), j = 2; h <= 1; j++) {
				results[i] = results[i] + h * pas / (1 - 1 / Math.PI * Math.acos(alpha[i] / h));
				h = alpha[i] + (j * pas);
			}

			results[i] = results[i] * 2 + alpha[i] * alpha[i];
		}

	}

	// crée la liste d'appareillement des spots les plus proches, avec distance
	private ArrayList<apparatedSpots> appDetectConstruction(List<Spot> spots, List<Spot> spots2) {
		ArrayList<apparatedSpots> liste_retour = new ArrayList<apparatedSpots>();

		if (spots.isEmpty()) {
		} else {
			int nbdeta = spots.size();
			int nbdetb = spots2.size();
			double x_a, y_a, x_b, y_b, z_a, z_b;
			for (int p2 = 0; p2 < nbdetb; p2++) {
				x_b = spots2.get(p2).mass_center.x;
				y_b = spots2.get(p2).mass_center.y;
				z_b = spots2.get(p2).mass_center.z;
				x_a = spots.get(0).mass_center.x;
				y_a = spots.get(0).mass_center.y;
				z_a = spots.get(0).mass_center.z;
				double min_dist = Math.sqrt(Math.pow(x_a - x_b, 2) + Math.pow(y_a - y_b, 2) + Math.pow(z_a - z_b, 2));
				int indice = 0;
				for (int p = 0; p < nbdeta; p++) {
					x_a = spots.get(p).mass_center.x;
					y_a = spots.get(p).mass_center.y;
					double temp = Math.sqrt(Math.pow(x_a - x_b, 2) + Math.pow(y_a - y_b, 2) + Math.pow(z_a - z_b, 2));
					if (temp < min_dist) {
						min_dist = temp;
						indice = p;
					}
				}
				apparatedSpots aS = new apparatedSpots(spots.get(indice), spots2.get(p2), min_dist);
				liste_retour.add(aS);
			}
		}
		return liste_retour;
	}

	// crée la sous liste des spots "colocalisés" par rapport au pourcentage
	// calculé statistiquement
	ArrayList<apparatedSpots> appDetectSelect(ArrayList<apparatedSpots> liste_app_detect, int ind_max) {
		ArrayList<apparatedSpots> liste_retour = new ArrayList<apparatedSpots>();
		// détermination de la distance "max" acceptée
		// création d'une liste de toutes les distance
		ArrayList<Double> distances = new ArrayList<Double>();
		for (apparatedSpots aS : liste_app_detect) {
			distances.add(aS.distance);
		}
		Collections.sort(distances);
		// détermination de l'indice max qui donne la distance max
		// int ind_max = (int)(percentage*distances.size());
		if (distances.isEmpty() || ind_max == 0) {
		} else {
			double distance_max = distances.get(ind_max - 1);
			// création de la sous liste des apparayedSpots dont la distance est
			// < distance_max
			for (apparatedSpots aS : liste_app_detect) {
				if (aS.distance <= distance_max) {
					liste_retour.add(aS);
				}
			}
		}
		return liste_retour;
	}

	// remplit un tableau de ROIs correspondant aux spots qui colocalisent à
	// partir de la sélection des spots colocalisés
	private void roiColoc3D(int t, ROI[] coloc1, ROI[] coloc2, List<Spot> detections1, List<Spot> detections2,
			ArrayList<apparatedSpots> spotsColoc) {
		int ind = 0;
		for (apparatedSpots aS : spotsColoc) {
			Spot s1 = aS.s1;
			Spot s2 = aS.s2;
			double x1 = s1.mass_center.x;
			double y1 = s1.mass_center.y;
			double z1 = s1.mass_center.z;
			double x2 = s2.mass_center.x;
			double y2 = s2.mass_center.y;
			double z2 = s2.mass_center.z;
			for (Spot spot1 : detections1) {
				double x_1 = spot1.mass_center.x;
				double y_1 = spot1.mass_center.y;
				double z_1 = spot1.mass_center.z;
				if ((x1 == x_1) & (y1 == y_1) & (z1 == z_1)) {
					int siz = spot1.point3DList.size();
					icy.type.point.Point3D[] li = new icy.type.point.Point3D[siz];
					int l = 0;
					for (plugins.nchenouard.spot.Point3D pt3 : spot1.point3DList) {
						icy.type.point.Point3D pt = new icy.type.point.Point3D.Double(pt3.x, pt3.y, pt3.z);
						li[l] = pt;
						l++;
					}
					BooleanMask3D bm = new BooleanMask3D(li);
					ROI3DArea roi1 = new ROI3DArea(bm);
					roi1.setT(t);
					coloc1[ind] = roi1;
					break;
				}
			}

			for (Spot spot2 : detections2) {
				double x_2 = spot2.mass_center.x;
				double y_2 = spot2.mass_center.y;
				double z_2 = spot2.mass_center.z;
				if ((x2 == x_2) & (y2 == y_2) & (z2 == z_2)) {
					int siz = spot2.point3DList.size();
					icy.type.point.Point3D[] li = new icy.type.point.Point3D[siz];
					int l = 0;
					for (plugins.nchenouard.spot.Point3D pt3 : spot2.point3DList) {
						icy.type.point.Point3D pt2 = new icy.type.point.Point3D.Double(pt3.x, pt3.y, pt3.z);
						li[l] = pt2;
						l++;
					}
					BooleanMask3D bm = new BooleanMask3D(li);
					ROI3DArea roi2 = new ROI3DArea(bm);
					roi2.setT(t);
					coloc2[ind] = roi2;
					break;
				}
			}
			ind++;
		}
	}

	private void roiColoc2D(int t, ROI[] coloc1, ROI[] coloc2, List<Spot> detections1, List<Spot> detections2,
			ArrayList<apparatedSpots> spotsColoc) {
		int ind = 0;
		for (apparatedSpots aS : spotsColoc) {
			Spot s1 = aS.s1;
			Spot s2 = aS.s2;
			double x1 = s1.mass_center.x;
			double y1 = s1.mass_center.y;
			double x2 = s2.mass_center.x;
			double y2 = s2.mass_center.y;
			for (Spot spot1 : detections1) {
				double x_1 = spot1.mass_center.x;
				double y_1 = spot1.mass_center.y;
				if ((x1 == x_1) & (y1 == y_1)) {
					int siz = spot1.point3DList.size();
					Point[] li = new Point[siz];
					int l = 0;
					for (plugins.nchenouard.spot.Point3D pt3 : spot1.point3DList) {
						Point pt = new Point((int) pt3.x, (int) pt3.y);
						li[l] = pt;
						l++;
					}
					BooleanMask2D bm = new BooleanMask2D(li);
					ROI2DArea roi1 = new ROI2DArea(bm);
					roi1.setT(t);
					coloc1[ind] = roi1;
					break;
				}
			}

			for (Spot spot2 : detections2) {
				double x_2 = spot2.mass_center.x;
				double y_2 = spot2.mass_center.y;
				if ((x2 == x_2) & (y2 == y_2)) {
					int siz = spot2.point3DList.size();
					Point[] li = new Point[siz];
					int l = 0;
					for (plugins.nchenouard.spot.Point3D pt3 : spot2.point3DList) {
						Point pt = new Point((int) pt3.x, (int) pt3.y);
						li[l] = pt;
						l++;
					}
					BooleanMask2D bm = new BooleanMask2D(li);
					ROI2DArea roi2 = new ROI2DArea(bm);
					roi2.setT(t);
					coloc2[ind] = roi2;
					break;
				}
			}
			ind++;

		}
	}

	public Point2D getMassCenter(ROI2D roi) {
		double x = 0, y = 0;
		long len = 0;

		final BooleanMask2D mask = roi.getBooleanMask(true);
		final boolean m[] = mask.mask;
		final int h = mask.bounds.height;
		final int w = mask.bounds.width;

		int off = 0;
		for (int j = 0; j < h; j++) {
			for (int i = 0; i < w; i++) {
				if (m[off++]) {
					x += i;
					y += j;
					len++;
				}
			}
		}

		final Point2D pos2d = roi.getPosition2D();
		return new Point2D.Double(pos2d.getX() + (x / len), pos2d.getY() + (y / len));
	}

	private void performAnalysis(EzVarSequence sequence_1, EzVarSequence sequence_2, EzVarInteger channel1,
			EzVarInteger channel2, EzVarSwimmingObject<DetectionResult> detections,
			EzVarSwimmingObject<DetectionResult> detections2, EzVarBoolean exportExcel) {

		// initialisation xls
		int row = 0;
		WritableWorkbook WW = null;
		WritableSheet WS = null;
		if (exportExcel.getValue()) {
			int page = 1;
			try {
				File f = exportExcelFile.getValue(true);
				if (!FileUtil.getFileExtension(f.getPath(), false).equalsIgnoreCase("xls"))
					f = new File(f.getPath() + ".xls");
				WW = XLSUtil.loadWorkbookForWrite(f);
			} catch (Exception e) {
				e.printStackTrace();
				return;
			}
			WS = XLSUtil.createNewPage(WW, "Page" + page);
			XLSUtil.setCellString(WS, 0, 0, "Date of XLS page:");
			row++;
			XLSUtil.setCellString(WS, 0, row, new Date().toString());
			row++;
		}

		Sequence sequence = sequence_1.getValue();
		int dim = 2;
		if (sequence.getSizeZ() > 1)
			dim = 3;

		// gestion des rois d'analyse en entrée
		if (list_roi.isEmpty()) {
			no_roi = true;
			for (int t = 0; t < sequence_1.getValue().getSizeT(); t++) {
				ROI roi = null;

				ROI2DRectangle r = new ROI2DRectangle(sequence.getBounds2D());
				for (int h = 0; h < sequence.getSizeZ(); h++) {
					r.setZ(h);
					r.setT(t);
					roi = r.getUnion(roi);
				}
				list_roi.add(roi);
			}
		} else {
			// on teste si les dimension en c sonts compatible
			double c1 = list_roi.get(0).getPosition5D().getC();
			for (ROI r : list_roi) {
				double c = r.getPosition5D().getC();
				if ((c != c1) && (c != (-1))) {
					new AnnounceFrame("ROI channels are incompatibles");
					return;
				}
			}

			// on teste si les dimensions en temps/z sont incompatibles pour une
			// union
			boolean one_z = false;
			boolean all_z = true;

			for (ROI r : list_roi) {
				if (r.getBounds5D().isInfiniteZ() == false) {
					all_z = false;
				}
				if (r.getBounds5D().isInfiniteZ()) {
					one_z = true;
				}
			}
			// gestion de l'exception
			if (one_z == true && all_z == false) {
				new AnnounceFrame("Incompatibility in Z dimensions between ROIs");
				return;
			}
		}

		// on construit emsuite la roi d'analyse au temps t=0
		ROI roi_t = ROI_t(dim, list_roi, 0, -1);
		int T;
		switch (method.getValue()) {
		case CORRELATION:
			if (exportExcel.getValue()) {

				XLSUtil.setCellString(WS, 0, row, "Time");
				XLSUtil.setCellString(WS, 1, row, "Pearson coefficient");
				XLSUtil.setCellString(WS, 2, row, "p value");
				XLSUtil.setCellString(WS, 3, row, "M1");
				XLSUtil.setCellString(WS, 4, row, "p value M1");
				XLSUtil.setCellString(WS, 5, row, "M2");
				XLSUtil.setCellString(WS, 6, row, "p value M2");
				XLSUtil.setCellString(WS, 7, row, "ICCS1");
				XLSUtil.setCellString(WS, 8, row, "ICCS2");
				row++;
			}

			// initialisation boucle en temps
			// temps total de la sequence
			if (sequence_1.getValue() != null) {
				T = sequence_1.getValue().getSizeT();
			} else {
				return;
			}

			// démarrage boucle en temps
			for (int t = 0; t < T; t += 1) {
				if (exportExcel.getValue()) {
					XLSUtil.setCellNumber(WS, 0, row, t);
				}
				roi_t = ROI_t(dim, list_roi, t, -1);

				double pearson[] = Correlation.pearson_TCL(sequence_1.getValue(), sequence_2.getValue(), t,
						channel1.getValue(), channel2.getValue(), roi_t);
				pearson_coeff.setValue(pearson[0]);
				pearson_coeff_label.setValue(StringUtil.toString(pearson[0], 2));

				pvalue_pearson.setValue(pearson[1]);
				pvalue_pearson_label.setValue(StringUtil.toString(pearson[1], 2));
				// gerer l'export excel
				if (exportExcel.getValue()) {
					XLSUtil.setCellNumber(WS, 1, row, pearson_coeff.getValue());
					XLSUtil.setCellNumber(WS, 2, row, pvalue_pearson.getValue());
				}
				ArrayList<ROI> list1 = sequence_1.getValue().getROIs();
				ArrayList<ROI> list2 = sequence_2.getValue().getROIs();
				// on enleve ensuite les ROIs "d'étude"
				ArrayList<ROI> list1_ = new ArrayList<ROI>();
				ArrayList<ROI> list2_ = new ArrayList<ROI>();
				for (ROI r : list1) {
					if (list_roi.contains(r)) {
					} else {
						list1_.add(r);
					}
				}
				for (ROI r : list2) {
					if (list_roi.contains(r)) {
					} else {
						list2_.add(r);
					}
				}
				// on trie les ROIs avec le bon canal
				ArrayList<ROI> list1__ = new ArrayList<ROI>();
				ArrayList<ROI> list2__ = new ArrayList<ROI>();
				for (ROI r : list1_) {
					if ((r.getPosition5D().getC() == channel1.getValue()) || (r.getBounds5D().isInfiniteC())) {
						list1__.add(r);
					}
				}
				for (ROI r : list2_) {
					if ((r.getPosition5D().getC() == channel2.getValue()) || (r.getBounds5D().isInfiniteC())) {
						list2__.add(r);
					}
				}

				//

				ROI roi1 = ROI_t(dim, list1__, t, channel1.getValue());
				ROI roi2 = ROI_t(dim, list2__, t, channel2.getValue());
				double manders[] = new double[4];
				if (roi1 == null || roi2 == null) {
					new AnnounceFrame("Incompatibility in Z dimensions between ROIs or no ROIs for Manders Analysis");
				} else {
					manders = Correlation.MandersCoeff(sequence_1.getValue(), sequence_2.getValue(), t,
							channel1.getValue(), channel2.getValue(), surface.getValue(), roi1, roi2, roi_t);
					M1.setValue(manders[0]);
					M1_label.setValue(StringUtil.toString(manders[0], 2));
					M2.setValue(manders[1]);
					M2_label.setValue(StringUtil.toString(manders[1], 2));

					pvalue_M1.setValue(manders[2]);
					pvalue_M1_label.setValue(StringUtil.toString(manders[2], 2));
					pvalue_M2.setValue(manders[3]);
					pvalue_M2_label.setValue(StringUtil.toString(manders[3], 2));

					if (exportExcel.getValue()) {
						XLSUtil.setCellNumber(WS, 3, row, M1.getValue());
						XLSUtil.setCellNumber(WS, 4, row, pvalue_M1.getValue());
						XLSUtil.setCellNumber(WS, 5, row, M2.getValue());
						XLSUtil.setCellNumber(WS, 6, row, pvalue_M2.getValue());
					}
				}

				double results_iccs[] = Correlation.ICCS_compute(sequence_1.getValue(), sequence_2.getValue(),
						channel1.getValue(), channel2.getValue(), t, roi_t);
				ICCS1.setValue(results_iccs[0]);
				ICCS1_label.setValue(StringUtil.toString(results_iccs[0], 2));

				ICCS2.setValue(results_iccs[1]);
				ICCS2_label.setValue(StringUtil.toString(results_iccs[1], 2));
				if (exportExcel.getValue()) {
					XLSUtil.setCellNumber(WS, 7, row, ICCS1.getValue());
					XLSUtil.setCellNumber(WS, 8, row, ICCS2.getValue());
				}
				row++;
			}

			break;

		case OBJECT:
			if (exportExcel.getValue()) {
				XLSUtil.setCellString(WS, 1, row, "Nb detections 1");
				XLSUtil.setCellString(WS, 2, row, "Nb detections 2");
				XLSUtil.setCellString(WS, 3, row, "r_1");
				XLSUtil.setCellString(WS, 4, row, "r_max");
				XLSUtil.setCellString(WS, 5, row, "Max. of the K function");
				{
					XLSUtil.setCellString(WS, 6, row, "p-value");
					XLSUtil.setCellString(WS, 7, row, "log_10(p-value)");
					XLSUtil.setCellString(WS, 8, row, "Percentage of detections 2 coloc. with detections 1");
					XLSUtil.setCellString(WS, 9, row, "Number of detections 2 coloc. with detections 1");
					XLSUtil.setCellString(WS, 10, row, "Distance of Colocalization");
					XLSUtil.setCellString(WS, 11, row, "Std. of Coloc. Dist.");
					XLSUtil.setCellString(WS, 13, row, "radius");
					XLSUtil.setCellString(WS, 14, row, "% of coloc");
				}
				row++;
			}

			// définition du ratio Z/X
			double ratio_zx = 1;
			if (sequence.getSizeZ() > 1) {
				ratio_zx = sequence.getPixelSizeZ() / sequence.getPixelSizeX();
			}
			// initialisation boucle en temps
			// temps total de la sequence

			DetectionResult detect_ = (DetectionResult) detections.getValue().getObject();
			sequence = detect_.getSequence();
			DetectionResult detect2_ = (DetectionResult) detections2.getValue().getObject();
			Sequence sequence2 = detect2_.getSequence();
			T = sequence.getSizeT();
			betaCorrection(max_radius.getValue() / 10, 100);
			/////////////////////////////////////////////
			//////////////////////////////////////////////////////
			if (hull) {
				// il faut reconstruire la liste de ROI
				list_roi.clear();
				for (ROI r : hull_roi) {
					// on efface les roi existantes sur la sequence
					sequence.removeROI(r);
					sequence2.removeROI(r);
				}
				hull_roi.clear();

				// boucle en temps
				for (int t = 0; t < T; t++) {
					// Vector<Spot> detection_0=detect_.getDetectionsAtT(t);
					Vector<Spot> detection2_00 = detect2_.getDetectionsAtT(t);
					List<Spot> detection2_0 = new ArrayList<Spot>();
					int index = 0;
					for (Spot s : detection2_00) {
						detection2_0.add(detection2_00.get(index));
						index++;
					}

					if (dim == 2) {
						List<Point2D> points = new ArrayList<Point2D>();
						/*
						 * for (Spot s:detection_0){
						 * plugins.nchenouard.spot.Point3D pt3 = s.mass_center;
						 * Point pt = new Point((int)pt3.x,(int)pt3.y);
						 * points.add(pt);}
						 */
						for (Spot s : detection2_0) {
							plugins.nchenouard.spot.Point3D pt3 = s.mass_center;
							Point pt = new Point((int) pt3.x, (int) pt3.y);
							points.add(pt);
						}

						List<Point2D> enveloppe = QuickHull2D.computeConvexEnvelope(points);
						ROI2DPolygon poly = new ROI2DPolygon(enveloppe);
						poly.setT(t);
						hull_roi.add(poly);
						list_roi.add(poly);
						sequence.addROI(poly);
						sequence2.addROI(poly);
					} else {
						Point3d[] points = new Point3d[detection2_0.size()];
						int ind = 0;
						/*
						 * for (Spot s:detection_0){ Point3d pt3 = new
						 * Point3d(s.mass_center.x, s.mass_center.y,
						 * s.mass_center.z); points[ind]=pt3;ind++;}
						 */
						for (Spot s : detection2_00) {
							Point3d pt3 = new Point3d(s.mass_center.x, s.mass_center.y, s.mass_center.z);
							points[ind] = pt3;
							ind++;
						}

						QuickHull3D hull = new QuickHull3D();
						hull.build(points);
						ROI3DPolygonalMesh roi = new ROI3DPolygonalMesh(hull);
						roi.setT(t);
						hull_roi.add(roi);
						list_roi.add(roi);
						sequence.addROI(roi);
						sequence2.addROI(roi);
					}
				}
			}
			///////////////////////////////////////////////////
			///////////////////////////////////////////////////
			// initialisation
			Vector<Spot> detection_00 = detect_.getDetectionsAtT(0);
			Vector<Spot> detection2_00 = detect2_.getDetectionsAtT(0);
			List<Spot> detection_0 = new ArrayList<Spot>();
			int index = 0;
			for (Spot s : detection_00) {
				detection_0.add(detection_00.get(index));
				index++;
			}
			List<Spot> detection2_0 = new ArrayList<Spot>();
			index = 0;
			for (Spot s : detection2_00) {
				detection2_0.add(detection2_00.get(index));
				index++;
			}
			// on construit emsuite la roi au temps t=0
			roi_t = ROI_t(dim, list_roi, 0, -1);

			List<Spot> detection = detectionsInRoi(detection_0, roi_t, 0);
			List<Spot> detection2 = detectionsInRoi(detection2_0, roi_t, 0);

			double volume = roi_t.getNumberOfPoints() * ratio_zx;
			if (roi_t.getBounds5D().isInfiniteZ()) {
				volume = volume * sequence.getSizeZ() * ratio_zx;
			}

			int nbdeta = detection.size();
			int nbdetb = detection2.size();
			// construction du distance_tab pour calculer la p_value de coloc
			double coeff_radius = 1;
			double step_min = max_radius.getValue() / 10;

			if (sequence.getSizeZ() > 1) {
				mindist = Math.max(Math.pow(coeff_radius * volume / (nbdeta * nbdetb), 0.333), step_min);
			} else {
				mindist = Math.max(Math.pow(coeff_radius * volume / (nbdeta * nbdetb), 0.5), step_min);
			}
			if (mindist > max_radius.getValue()) {
				mindist = max_radius.getValue();
				new AnnounceFrame("Number of spots should be insufficient for the statistical analysis");
			}
			distance_fit.clear();
			distance_fit.add((double) 0);
			distance_fit.add(mindist);
			double temp = mindist;
			if (sequence.getSizeZ() > 1) {
				while (Math.max(Math.pow(Math.pow(temp, 3) + coeff_radius * volume / (nbdeta * nbdetb), 0.333),
						temp + step_min) < max_radius.getValue()) {
					temp = Math.max(Math.pow(Math.pow(temp, 3) + coeff_radius * volume / (nbdeta * nbdetb), 0.333),
							temp + step_min);
					distance_fit.add(temp);
				}
			} else {
				while (Math.max(Math.sqrt(Math.pow(temp, 2) + coeff_radius * volume / (nbdeta * nbdetb)),
						temp + step_min) < max_radius.getValue()) {
					temp = Math.max(Math.sqrt(Math.pow(temp, 2) + coeff_radius * volume / (nbdeta * nbdetb)),
							temp + step_min);
					distance_fit.add(temp);
				}
			}
			int N_fit = distance_fit.size();
			////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////
			/*
			 * double step = max_radius.getValue() / 10; distance_fit.clear();
			 * distance_fit.add((double) 0); distance_fit.add(step); double
			 * temp_fit = step; while (temp_fit + step < max_radius.getValue())
			 * { temp_fit += step; distance_fit.add(temp_fit); } int N_fit =
			 * distance_fit.size(); //pour le plot.. double[]dist_fit = new
			 * double[N_fit]; for (int g=0;g<N_fit;g++){
			 * dist_fit[g]=distance_fit.get(g); }
			 */
			////////////////////////////////////////
			if (export_colocalized_rois.getValue()) {
				if (coloc_roi_1 != null) {
					for (ROI roi : coloc_roi_1) {
						sequence.removeROI(roi);
					}
				}
				if (coloc_roi_2 != null) {
					for (ROI roi : coloc_roi_2) {
						sequence2.removeROI(roi);
					}
				}
			}

			// stockage de l'histogram des distances de coloc
			double[][] histo_dist = new double[T][N_fit - 1];
			// démarrage boucle en temps
			for (int t = 0; t < T; t += 1) {
				if (exportExcel.getValue()) {
					XLSUtil.setCellNumber(WS, 0, row, t);
				}

				// initialisation
				detection_0 = detect_.getDetectionsAtT(t);
				detection2_0 = detect2_.getDetectionsAtT(t);
				roi_t = ROI_t(dim, list_roi, t, -1);
				detection = detectionsInRoi(detection_0, roi_t, t);
				detection2 = detectionsInRoi(detection2_0, roi_t, t);

				// calcul des parametres globaux: aire totale des rois et nb de
				// detections
				// faut-refaire une fonction calcul de volume
				if (roi_t == null) {
					volume = 0;
				} else {
					volume = roi_t.getNumberOfPoints() * ratio_zx;
					if (roi_t.getBounds5D().isInfiniteZ()) {
						volume = volume * sequence.getSizeZ() * ratio_zx;
					}
				}

				nbdeta = detection.size();
				nbdetb = detection2.size();
				// initialisation
				double[] delta_K = new double[N_fit - 1];
				double[][] Ktemp = new double[N_fit - 1][3];
				double[] coeffs = new double[3];
				double[] var = new double[N_fit];
				double[] retour = new double[5];

				if (manual.getValue()) {
					// calcul du pourcentage de molecules avec distancea< rmax
					double alpha_manual = manual_computation(detection, detection2, max_radius);
					alpha_fit.setValue(alpha_manual);
					alpha_fit_label.setValue(StringUtil.toString(alpha_manual, 2));

					number_fit.setValue(Math.ceil(alpha_manual * nbdetb));
					number_fit_label.setValue(StringUtil.toString(Math.ceil(alpha_manual * nbdetb), 2));
				} else {
					////////////////////////////////////////////////////////////
					////////////////////////////////////////////////////////////
					////////////////////////////////////////////////////////////
					////////////////////////////////////////////////////////////
					////////////////////////////////////////////////////////////
					////////////////////////////////////////////////////////////
					Chronometer chrono = new Chronometer("chrono");

					if (sequence.getSizeZ() > 1) {
						Ktemp = Ripley3D.correlation(sequence.getSizeZ(), roi_t, detection, detection2, volume, nbdeta,
								nbdetb, N_fit, distance_fit, ratio_zx);
					} else {
						Ktemp = Ripley2D.correlation(roi_t, detection, detection2, volume, nbdeta, nbdetb, N_fit,
								distance_fit);
					}
					chrono.displayMs();

					for (int k = 0; k < N_fit - 1; k++) {
						delta_K[k] = Ktemp[k][0];
					}

					double[] Kfit = new double[N_fit];
					Kfit[0] = 0;
					for (int i = 0; i < N_fit - 1; i++) {
						Kfit[i + 1] = Kfit[i] + delta_K[i];
					}

					if (sequence.getSizeZ() > 1) {
						var = Ripley3D.variance_theo(sequence.getSizeZ(), roi_t, detection, detection2, volume, nbdeta,
								nbdetb, N_fit, distance_fit, N_h, ratio_zx);
					} else
					// il faut calculer la roi au temps t//construire une array
					// list
					{
						var = Ripley2D.variance_theo(roi_t, detection, detection2, volume, nbdeta, nbdetb, N_fit,
								distance_fit, N_h, results);
					}
					K.clear();
					K.add(0, 0.0);
					// calcul du sup de la fonction K normalisée
					double s = -10000;
					for (int h = 1; h < N_fit; h++) {
						if (sequence.getBounds5D().getSizeZ() > 1) {
							Kfit[h] = Kfit[h] - (double) (4 / 3) * Math.PI * Math.pow(distance_fit.get(h), 3);
						} else {
							Kfit[h] = Kfit[h] - Math.PI * Math.pow(distance_fit.get(h), 2);
						}
						K.add(h, Kfit[h] / Math.sqrt(var[h]));
						if (K.get(h) > s)
							s = K.get(h);
					}

					sup_K.setValue(s);
					sup_K_label.setValue(StringUtil.toString(s, 2));
					double normcdf = 0.5 * (1 + ErrorFunction.erf(s / Math.sqrt(2)));
					p_value_K.setValue(1 - normcdf);
					if (s < 4) {
						p_value_K_label.setValue(StringUtil.toString((double) p_value_K.getValue(), 4));
					} else // si s>4
					{
						double x = s / Math.sqrt(2);
						double y = 1 / (2 * Math.sqrt(Math.PI) * x);
						double exponent = (Math.log(y) - Math.pow(x, 2)) / Math.log(10);
						p_value_K_label.setValue("10^ " + StringUtil.toString(exponent, 2));
					}
					// fit sur K_fit moyen

					double[] start_est = new double[3];
					start_est[0] = 0.1D; // initial estimate of alpha
					start_est[1] = 1.01D; // initial estimate of mu
					start_est[2] = 0.3D;// sigma_ini.getValue(); // initial
										// estimate of sigma

					// voir pour entrer u
					// double sigma_fit=0.3;
					// fit_data.main(distance_fit, Kfit, area, nbdeta,
					// start_est);
					if (sequence.getBounds5D().getSizeZ() > 1) {
						fit_data.main_manual(distance_fit, Kfit, volume, nbdeta);
						coeffs = fit_data.coeffs;
					} else {
						fit_data.main_manual(distance_fit, Kfit, volume, nbdeta);
						coeffs = fit_data.coeffs;
					}

					alpha_fit.setValue(coeffs[0]);
					alpha_fit_label.setValue(StringUtil.toString(coeffs[0], 2));

					mu_fit.setValue(coeffs[1]);
					mu_fit_label.setValue(StringUtil.toString(coeffs[1], 2));

					number_fit.setValue(Math.ceil(alpha_fit.getValue() * nbdetb));
					number_fit_label.setValue(StringUtil.toString(number_fit.getValue(), 2));

					sigma_fit.setValue(coeffs[2]);
					sigma_fit_label.setValue(StringUtil.toString(coeffs[2], 2));
				}
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

				if (export_colocalized_rois.getValue()) {

					// Une fois que l'on connait le pourcentage de coloc, on
					// peut construire une fonction qui extrait les ROIs
					// correspondanst aux spots coloc
					// 1 fonction pour déterminer la liste des apparatedSpots
					ArrayList<apparatedSpots> liste_app_detect = new ArrayList<apparatedSpots>();
					liste_app_detect = appDetectConstruction(detection, detection2);
					// 1 fonction pour construire la listes des detections 1 et
					// 2 qui colocalisent à partir de la liste des
					// apparatedDetection et du pourcentage de coloc
					int ind_max = (int) (Math.ceil(alpha_fit.getValue() * nbdetb));
					ArrayList<apparatedSpots> liste_short = new ArrayList<apparatedSpots>();
					liste_short = appDetectSelect(liste_app_detect, ind_max);
					// 1 fonction (utilisée 2 fois) pour transformer les 2
					// listes de detections en listes de ROIs
					// création des 2 vecteurs de ROIs
					int si = liste_short.size();
					coloc_roi_1 = new ROI[si];
					coloc_roi_2 = new ROI[si];

					if (sequence.getSizeZ() > 1) {
						roiColoc3D(t, coloc_roi_1, coloc_roi_2, detection, detection2, liste_short);
					} else {
						roiColoc2D(t, coloc_roi_1, coloc_roi_2, detection, detection2, liste_short);
					}
					for (ROI roi : coloc_roi_1) {
						sequence.addROI(roi);
					}
					for (ROI roi : coloc_roi_2) {
						sequence2.addROI(roi);
					}

				}

				if (exportExcel.getValue()) {
					XLSUtil.setCellNumber(WS, 1, row, nbdeta);
					XLSUtil.setCellNumber(WS, 2, row, nbdetb);
					XLSUtil.setCellNumber(WS, 3, row, mindist);
					XLSUtil.setCellNumber(WS, 4, row, max_radius.getValue());
					XLSUtil.setCellNumber(WS, 5, row, sup_K.getValue());
					XLSUtil.setCellNumber(WS, 6, row, p_value_K.getValue());
					if (sup_K.getValue() < 4) {
						XLSUtil.setCellNumber(WS, 7, row, Math.log10(p_value_K.getValue()));
					} else {
						double x = sup_K.getValue() / Math.sqrt(2);
						double y = 1 / (2 * Math.sqrt(Math.PI) * x);
						double exponent = (Math.log(y) - Math.pow(x, 2)) / Math.log(10);
						XLSUtil.setCellNumber(WS, 7, row, exponent);
					}
					XLSUtil.setCellNumber(WS, 8, row, alpha_fit.getValue());
					XLSUtil.setCellNumber(WS, 9, row, number_fit.getValue());
					XLSUtil.setCellNumber(WS, 10, row, mu_fit.getValue());
					XLSUtil.setCellNumber(WS, 11, row, sigma_fit.getValue());
					for (int i = 0; i < N_fit - 1; i++) {
						XLSUtil.setCellNumber(WS, 13, row, distance_fit.get(i + 1));
						XLSUtil.setCellNumber(WS, 14, row, histo_dist[t][i]);
						row++;
					}
					row++;
				}
				if (graph.getValue()) {

					K_fit.clear();
					K_fit.add(0, 0.0);
					for (int p = 1; p < N_fit; p++) {
						K_fit.add(p,
								(alpha_fit.getValue() * volume * 0.5 / nbdeta
										* (Stat.erf((distance_fit.get(p) - mu_fit.getValue())
												/ (sigma_fit.getValue() * Math.sqrt(2)))
												- Stat.erf((-distance_fit.get(p) - mu_fit.getValue())
														/ (sigma_fit.getValue() * Math.sqrt(2)))))
										/ (Math.sqrt(var[p])));
					}
					plotGraph(distance_fit, K, K_fit);
				}
			}
			if (no_roi)
				list_roi.clear();
		}
		///////////////////////////////
		if (exportExcel.getValue()) {
			try {
				XLSUtil.saveAndClose(WW);
			} catch (WriteException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

	}

	private double manual_computation(List<Spot> detection, List<Spot> detection2, EzVarDouble max_radius2) {
		int compteur = 1;
		ArrayList<Spot> done = new ArrayList<Spot>();
		for (Spot s1 : detection) {
			for (Spot s2 : detection2) {
				double distance_12 = Math.pow(s1.mass_center.x - s2.mass_center.x, 2)
						+ Math.pow(s1.mass_center.y - s2.mass_center.y, 2);
				if (distance_12 < Math.pow(max_radius2.getValue(), 2)) {
					if (done.contains(s2)) {
					} else {
						done.add(s2);
						compteur += 1;
					}
				}
			}
		}
		double temp = Math.pow(detection2.size(), -1);
		double alpha = compteur * temp;
		return alpha;
	}

	@Override
	public void clean() {
		// TODO Auto-generated method stub

	}

	@Override
	protected void execute() {
		// double pas = 0.05;
		// int N = 5;

		pearson_coeff_label.setValue("  N/A  ");
		M1_label.setValue("  N/A  ");
		M2_label.setValue("  N/A  ");
		pvalue_pearson_label.setValue("  N/A  ");
		pvalue_M2_label.setValue("  N/A  ");
		pvalue_M1_label.setValue("  N/A  ");
		ICCS1_label.setValue("  N/A  ");
		ICCS2_label.setValue("  N/A  ");
		sup_K_label.setValue("  N/A  ");
		alpha_fit_label.setValue("  N/A  ");
		number_fit_label.setValue("  N/A  ");
		mu_fit_label.setValue("  N/A  ");
		sigma_fit_label.setValue("  N/A  ");
		p_value_K_label.setValue("  N/A  ");
		switch (method.getValue()) {
		case CORRELATION:
			if (sequence1.getValue() == null || sequence2.getValue() == null) {
				new AnnounceFrame("Please first select sequences");
				return;
			}
			if (sequence1.getValue().getSizeX() != sequence2.getValue().getSizeX()
					|| sequence1.getValue().getSizeY() != sequence2.getValue().getSizeY()
					|| sequence1.getValue().getSizeZ() != sequence2.getValue().getSizeZ()
					|| sequence1.getValue().getSizeT() != sequence2.getValue().getSizeT()) {
				new AnnounceFrame("Sequences must have the same (X,Y,Z,T) dimensions");
				return;
			}
			if (channel1.getValue() > sequence1.getValue().getSizeC() - 1
					|| channel2.getValue() > sequence2.getValue().getSizeC() - 1) {
				new AnnounceFrame("Channels do not exist");
				return;
			}
			performAnalysis(sequence1, sequence2, channel1, channel2, detections1, detections2, exportExcel);
			break;
		case OBJECT:
			if (max_radius.getValue() <= 0) {
				new AnnounceFrame("Max. radius must be > 0");
				return;
			}
			if (detections1.getValue() == null || detections2.getValue() == null) {
				new AnnounceFrame("Please first select detections set");
				return;
			}
			if (sequence1.getValue() == null || sequence2.getValue() == null) {
				new AnnounceFrame("There is no sequence associated to detections");
				return;
			} else {
				performAnalysis(sequence1, sequence2, channel1, channel2, detections1, detections2, exportExcel);
			}
			break;
		}
	}

	@Override
	protected void initialize() {

		super.addEzComponent(method);

		method.addVisibilityTriggerTo(hull_b, ColocMethod.OBJECT);
		super.addEzComponent(hull_b);

		super.addEzComponent(setTrackSetA);
		super.addEzComponent(reset);

		EzGroup sequences = new EzGroup("Sequences for correlation analysis", sequence1, channel1, sequence2, channel2);
		method.addVisibilityTriggerTo(sequences, ColocMethod.CORRELATION);
		// method.addVisibilityTriggerTo(rois, ColocMethod.CORRELATION);
		super.addEzComponent(sequences);

		// super.addEzComponent(rois);

		EzGroup detections = new EzGroup("Spot detections for object-based analysis", detections1, detections2);
		method.addVisibilityTriggerTo(detections, ColocMethod.OBJECT);
		addEzComponent(detections);

		// min_radius = new EzVarDouble("Minimal distance for the analysis (in
		// pxs)", 0, 0, 100, 1);
		max_radius = new EzVarDouble("Maximal distance for the analysis (in pxs)", 5, 0, 100, 1);
		graph = new EzVarBoolean("Plot K function graph to check the fit", false);
		manual = new EzVarBoolean("Manual estimation", false);

		// sigma_ini = new EzVarDouble("Initialization of coloc distance
		// std.",0.3, 0, 3, 0.01);

		// EzGroup K_param = new
		// EzGroup("Parameters",max_radius,sigma_ini,graph);
		EzGroup K_param = new EzGroup("Parameters", max_radius, graph, manual);

		method.addVisibilityTriggerTo(K_param, ColocMethod.OBJECT);
		addEzComponent(K_param);

		addEzComponent(new EzGroup("Export", export_colocalized_rois, exportExcel, exportExcelFile));
		method.addVisibilityTriggerTo(export_colocalized_rois, ColocMethod.OBJECT);
		exportExcel.addVisibilityTriggerTo(exportExcelFile, true);

		/////////////////////////////////////////////////////////////
		addComponent(new JSeparator(JSeparator.VERTICAL));
		////////////////////////////////////////////////////////////////

		EzGroup pears = new EzGroup("Pearson Analysis: Results", pearson_coeff_label, pvalue_pearson_label);
		method.addVisibilityTriggerTo(pears, ColocMethod.CORRELATION);
		super.addEzComponent(pears);

		EzGroup manders = new EzGroup("Manders Analysis: Results", M1_label, pvalue_M1_label, M2_label,
				pvalue_M2_label);
		method.addVisibilityTriggerTo(manders, ColocMethod.CORRELATION);
		super.addEzComponent(manders);

		// EzGroup icq = new EzGroup("ICQ Analysis:
		// Results",ICQ_label,pvalue_ICQ_label);
		// method.addVisibilityTriggerTo(icq, ColocMethod.CORRELATION);
		// super.addEzComponent(icq);

		EzGroup iccs = new EzGroup("Cross-Correlation Analysis: Results", ICCS1_label, ICCS2_label);
		method.addVisibilityTriggerTo(iccs, ColocMethod.CORRELATION);
		super.addEzComponent(iccs);

		sup_K = new EzVarDouble("Max. of the K function", 0, 0, 1, 0.01);
		p_value_K = new EzVarDouble("p-value", 0, 0, 1, 0.01);
		alpha_fit = new EzVarDouble("Percentage of detections 2 colocalizing with detections 1 (fit)", 0.1, 0, 1, 0.01);
		number_fit = new EzVarDouble("Number of detections 2 colocalizing with detections 1 (fit)", 1, 0, 1, 0.1);
		mu_fit = new EzVarDouble("Distance of coloc. (in pixels) (fit)", 1, 0, 10, 0.1);
		sigma_fit = new EzVarDouble("Std. of coloc. distance", 0.3, 0, 3, 0.01);
		EzGroup K_results = new EzGroup("Object-based analysis: Results", sup_K_label, p_value_K_label, alpha_fit_label,
				number_fit_label, mu_fit_label, sigma_fit_label);

		method.addVisibilityTriggerTo(K_results, ColocMethod.OBJECT);
		addEzComponent(K_results);
		////////////////////////////////////////////////////////////
		addComponent(new JSeparator(JSeparator.VERTICAL));
		////////////////////////////////////////////////////////////////

	}

	private void plotGraph(ArrayList<Double> dist_tab, ArrayList<Double> measurementList,
			ArrayList<Double> fittedList) {

		final JFreeChart chart;
		XYSeriesCollection xyDataset = new XYSeriesCollection();

		chart = ChartFactory.createXYLineChart("Ripley's K function step", "pixels", "K", xyDataset,
				PlotOrientation.VERTICAL, true, false, false);

		XYSeries seriesXY = new XYSeries("Y");
		xyDataset.addSeries(seriesXY);
		XYSeries seriesXYfit = new XYSeries("Y fit");
		xyDataset.addSeries(seriesXYfit);

		for (int t = 0; t < measurementList.size(); t++) {
			seriesXY.add((double) dist_tab.get(t), measurementList.get(t));
			seriesXYfit.add((double) dist_tab.get(t), fittedList.get(t));
		}

		if (this.chartpanel == null) {
			this.chartpanel = new ChartPanel(chart);
			addComponent(this.chartpanel);
			getUI().repack(true);
		} else {
			this.chartpanel.setChart(chart);
		}

		/*
		 * ThreadUtil.invokeLater(new Runnable() {
		 * 
		 * @Override public void run() { IcyFrame graphFrame = new
		 * IcyFrame("graph", true, true, true, true);
		 * graphFrame.setContentPane(new ChartPanel(chart, 500, 200, 500, 200,
		 * 500, 500, false, false, true, true, true, true));
		 * graphFrame.setVisible(true); graphFrame.pack();
		 * graphFrame.addToMainDesktopPane(); graphFrame.center(); } });
		 */

	}

	@Override
	public void sequenceChanged(SequenceEvent sequenceEvent) {
		// TODO Auto-generated method stub

	}

	@Override
	public void sequenceClosed(Sequence sequence) {
		// TODO Auto-generated method stub

	}

	@Override
	public void actionPerformed(ActionEvent e) {
		if (e.getActionCommand() == setTrackSetA.name) {
			list_roi.clear();
			for (ROI roi : Icy.getMainInterface().getROIs()) {
				if (roi.isSelected()) {
					list_roi.add(roi);
				}
			}
			setTrackSetA.setText("You have selected " + list_roi.size() + " ROIs");
			hull_b.setEnabled(false);
			hull = false;
		}
		if (e.getActionCommand() == hull_b.name) {
			list_roi.clear();
			hull = true;
			setTrackSetA.setEnabled(false);
			hull_b.setText("ROI = Detections 2 Convex hull");

		}

		if (e.getActionCommand() == reset.name) {
			list_roi.clear();
			hull = false;
			setTrackSetA.setText(setTrackSetA.name);
			setTrackSetA.setEnabled(true);
			hull_b.setText(hull_b.name);
			hull_b.setEnabled(true);
			for (ROI r : hull_roi) {
				for (Sequence seq : Icy.getMainInterface().getSequences()) {
					seq.removeROI(r);
				}
			}
			hull_roi.clear();
		}
	}

	private List<Spot> detectionsInRoi(List<Spot> detection, ROI roi, int t) {
		// SEQUENCE A remplacer par w,h

		// create a hashMap with the detections binded to ROI

		List<Spot> ROIDetection = new ArrayList<Spot>();
		if (roi == null) {
			return ROIDetection;
		}

		ArrayList<Spot> TestList = new ArrayList<Spot>();

		// fill hashMap
		for (Spot spot : detection) {
			boolean b = false;
			for (Spot spot2 : TestList) {
				if ((spot2.mass_center.x == spot.mass_center.x) && (spot2.mass_center.y == spot.mass_center.y)
						&& (spot2.mass_center.z == spot.mass_center.z)) {
					b = true;
				}
			}
			if (b == false)// si le spot n'a pas encore été pris en compte
			{
				TestList.add(spot);
				// if (roi.contains(spot.mass_center.x, spot.mass_center.y,
				// spot.mass_center.z, t, roi.getPosition5D().getC()))
				if (roi.contains(spot.mass_center.x, spot.mass_center.y, spot.mass_center.z, roi.getPosition5D().getT(),
						roi.getPosition5D().getC())) {
					ROIDetection.add(spot);
				}
			}
		}

		if (detection.size() == 0) {
			new AnnounceFrame("There is no detection associated with the ROI(s)");

		}

		return (ROIDetection);
	}

	private ROI ROI_t(int dim, ArrayList<ROI> list_roi_0, int t, int c) {
		ArrayList<ROI> list_roi = new ArrayList<ROI>();
		for (ROI ro : list_roi_0) {
			ROI ro2 = ro.getCopy();
			list_roi.add(ro2);
		}
		// check compatibility en z si dim =3
		if (dim == 3) {
			boolean one_z = false;
			boolean all_z = true;
			for (ROI r0 : list_roi) {
				if (r0.getBounds5D().isInfiniteZ() == false) {
					all_z = false;
				}
				if (r0.getBounds5D().isInfiniteZ()) {
					one_z = true;
				}
			}
			// gestion de l'exception
			if (one_z == true && all_z == false) {
				return null;
			}
		}

		ROI roi_t = null;

		for (ROI r : list_roi) {
			Point5D pt0 = r.getPosition5D();
			if (dim == 2) {
				pt0.setZ(-1);
			}
			if ((r.getBounds5D().isInfiniteC()) || (c == -1)) {
				pt0.setC(c);
			}
			if ((r.getBounds5D().isInfiniteT()) || (t == -1)) {
				pt0.setT(t);
			}
			r.setPosition5D(pt0);
			if ((pt0.getT() == t) && (pt0.getC() == c)) {
				roi_t = r.getUnion(roi_t);
			}
		}
		return roi_t;
	}

}
