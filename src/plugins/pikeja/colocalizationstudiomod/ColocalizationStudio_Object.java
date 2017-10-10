package plugins.pikeja.colocalizationstudiomod;

import java.awt.Point;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.List;
import java.util.Vector;

import org.apache.poi.hssf.usermodel.HSSFPalette;
import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.ss.usermodel.CellStyle;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.ss.util.WorkbookUtil;

import javassist.expr.Instanceof;

import icy.file.FileUtil;
import icy.gui.frame.progress.AnnounceFrame;
import icy.plugin.abstract_.Plugin;
import icy.plugin.interface_.PluginBundled;
import icy.roi.BooleanMask2D;
import icy.roi.BooleanMask3D;
import icy.roi.ROI;
import icy.roi.ROI2D;
import plugins.kernel.roi.roi2d.ROI2DPoint;
import icy.sequence.Sequence;
import icy.type.point.Point3D;
import icy.type.point.Point5D;
import icy.util.XLSUtil;
import jxl.write.WritableSheet;
import jxl.write.WritableWorkbook;
import jxl.write.WriteException;
import plugins.adufour.blocks.lang.Block;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.vars.lang.VarBoolean;
import plugins.adufour.vars.lang.VarDouble;

import plugins.adufour.vars.lang.VarFile;
import plugins.adufour.vars.lang.VarROIArray;
import plugins.adufour.vars.lang.VarSequence;
import plugins.adufour.vars.lang.VarWorkbook;
import plugins.kernel.roi.roi2d.ROI2DArea;
import plugins.kernel.roi.roi2d.ROI2DRectangle;
import plugins.kernel.roi.roi3d.ROI3DArea;
import plugins.nchenouard.spot.DetectionResult;
import plugins.nchenouard.spot.Spot;

// Colocalisation with Ripley function K
// Significant 

public class ColocalizationStudio_Object extends Plugin implements Block, PluginBundled {

	// EzVarSwimmingObject<DetectionResult> detections = new
	// EzVarSwimmingObject<DetectionResult>("Detections");
	VarROIArray detections1 = new VarROIArray("List of detections 1 (ROIS)");
	VarROIArray detections2 = new VarROIArray("List of detections 2 (ROIS)");

	// ROIs (spots) qui colocalisent en sortie
	VarBoolean export_colocalized_rois = new VarBoolean("Export colocalized detections (ROIs)", false);
	VarROIArray detections1Coloc = new VarROIArray("List of colocalized detections 1 (ROIS)");
	VarROIArray detections2Coloc = new VarROIArray("List of colocalized detections 2 (ROIS)");

	VarBoolean external_rois = new VarBoolean("Choose ROIs", false);
	VarROIArray ROIs = new VarROIArray("List of ROIs");

	VarSequence input_sequence1 = new VarSequence("Input sequence 1", null);
	VarSequence input_sequence2 = new VarSequence("Input sequence 2", null);

	VarDouble maxdistance_in = new VarDouble("Max. Radius (in pixels)", 5);
	// VarDouble min_step = new VarDouble("Min. Step (in pixels)", 1);

	// Estimation des parametres (pourcentage coloc et distance coloc)
	VarDouble number_fit = new VarDouble("Number of detections 2 colocalizing with detections 1", 0.1);
	VarDouble alpha_fit = new VarDouble("Percentage of detections 2 colocalizing with detections 1", 0.1);
	VarDouble mu_fit = new VarDouble("Distance of coloc. (in pixels)", 0.1);
	VarDouble sigma_fit = new VarDouble("Std. Dev. of Coloc. Distance (in pixels)", 0.1);

	// VarString pvalue = new VarString("p value", "");
	VarDouble pvalue = new VarDouble("p-value", 0.1);

	// VarBoolean exportExcel = new VarBoolean("Export results in excel",false);
	// VarFile exportExcelFile = new VarFile("Excel File", null);
	VarWorkbook book = new VarWorkbook("Workbook", (Workbook) null);

	double[] results = null;
	int N_h;
	// pour le test stat
	double maxdist;
	double mindist;
	ArrayList<Double> distance_fit = new ArrayList<Double>();
	ArrayList<ROI> list_roi = new ArrayList<ROI>();
	private ArrayList<Double> K = new ArrayList<Double>();

	@Override
	public void declareInput(VarList inputMap) {

		inputMap.add("Detections 1 (ROIs)", detections1);
		inputMap.add("Detections 2 (ROIs)", detections2);
		inputMap.add("Input Sequence 1", input_sequence1);
		inputMap.add("Input Sequence 2", input_sequence2);
		inputMap.add("Choose ROIs", external_rois);
		inputMap.add("ROIs", ROIs);
		inputMap.add("Max. Radius (in x pixels)", maxdistance_in);
		inputMap.add("Export colocalized spots (ROIs)", export_colocalized_rois);
		// inputMap.add( "Excel Export", exportExcel );
		// inputMap.add( "Excel File", exportExcelFile );
	}

	@Override
	public void declareOutput(VarList outputMap) {
		outputMap.add("percentage of colocalized detections 2", alpha_fit);
		outputMap.add("number of colocalized detections 2", number_fit);
		outputMap.add("distance of colocalized detections 2", mu_fit);
		outputMap.add("p value", pvalue);
		outputMap.add("Colocalized Detections 1", detections1Coloc);
		outputMap.add("Colocalized Detections 2", detections2Coloc);
		outputMap.add("Workbook", book);
	}
	// outputMap.add("max distance",maxdist);

	// N: nb entre dmin et dmax avec pas donne
	@Override
	public void run() {
		apparatedSpots.coloc1.clear();
		apparatedSpots.coloc2.clear();
		if (maxdistance_in.getValue() <= 0) {
			new AnnounceFrame("Max. radius must be > 0");
			return;
		}
		if (detections1 == null || detections2 == null) {
			new AnnounceFrame("Please first select a detection set");
			return;
		} else {
			DetectionResult detect_ = new DetectionResult();
			DetectionResult detect2_ = new DetectionResult();

			for (ROI roi : detections1.getValue()) {
				Point3D p = massCenters.getMassCenter(roi);
				Spot spot = new Spot(p.getX(), p.getY(), p.getZ());
				if (roi.getBounds5D().isInfiniteT()) {
					detect_.addDetection(-1, spot);
				} else {
					detect_.addDetection((int) roi.getPosition5D().getT(), spot);
				}
			}

			for (ROI roi : detections2.getValue()) {
				Point3D p = massCenters.getMassCenter(roi);
				Spot spot = new Spot(p.getX(), p.getY(), p.getZ());
				if (roi.getBounds5D().isInfiniteT()) {
					detect2_.addDetection(-1, spot);
				} else {
					detect2_.addDetection((int) roi.getPosition5D().getT(), spot);
				}
			}

			maxdist = maxdistance_in.getValue();
			if (input_sequence1.getValue() == null || input_sequence2.getValue() == null) {
				new AnnounceFrame("There is no open sequence");
				return;
			}

			else {
				detect_.setSequence(input_sequence1.getValue());
				detect2_.setSequence(input_sequence2.getValue());

				// créaton du workbook
				if (book.getValue() == null) {
					book.setValue(new HSSFWorkbook());
				}
				book.getValue().setMissingCellPolicy(Row.CREATE_NULL_AS_BLANK);

				performAnalysis(detect_, detect2_);
			}
		}
	}

	private void performAnalysis(DetectionResult detect_, DetectionResult detect2_) {

		/*
		 * // initialisation xls int row = 0; WritableWorkbook WW = null;
		 * WritableSheet WS = null; if (exportExcel.getValue()) { int page = 1;
		 * try { File f = exportExcelFile.getValue(true); if
		 * (!FileUtil.getFileExtension(f.getPath(), false)
		 * .equalsIgnoreCase("xls")) f = new File(f.getPath() + ".xls"); WW =
		 * XLSUtil.loadWorkbookForWrite(f); } catch (Exception e) {
		 * e.printStackTrace(); return;} WS = XLSUtil.createNewPage(WW, "Page" +
		 * page); XLSUtil.setCellString(WS, 0, 0, "Date of XLS page:"); row++;
		 * XLSUtil.setCellString(WS, 0, row, new Date().toString()); row++;}
		 */

		// initialisation du workbook
		Workbook wb = book.getValue();
		// create the sheet
		String sheetName = "Coloc. Object Analysis";
		Sheet sheet = wb.getSheet(sheetName);
		if (sheet == null)
			sheet = wb.createSheet(sheetName);

		Sequence sequence = input_sequence1.getValue();
		int dim = 2;
		if (sequence.getSizeZ() > 1)
			dim = 3;

		// gestion des rois d'analyse en entrée

		if (external_rois.getValue() == false) {
			for (int t = 0; t < input_sequence1.getValue().getSizeT(); t++) {
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
			list_roi.clear();
			ROI[] roiArray = ROIs.getValue();// TODO Ou sont les ROIs???
			if (roiArray.length == 0) {
				for (int t = 0; t < input_sequence1.getValue().getSizeT(); t++) {
					ROI roi = null;

					ROI2DRectangle r = new ROI2DRectangle(sequence.getBounds2D());
					for (int h = 0; h < sequence.getSizeZ(); h++) {
						r.setZ(h);
						r.setT(t);
						roi = r.getUnion(roi);
					}
					list_roi.add(roi);
				}
			} else
				for (ROI roi : roiArray) {
					list_roi.add(roi);
				}

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
		/*
		 * if (exportExcel.getValue()) { XLSUtil.setCellString(WS,
		 * 1,row,"Nb detections 1"); XLSUtil.setCellString(WS,
		 * 2,row,"Nb detections 2"); XLSUtil.setCellString(WS, 3,row,"r_1");
		 * XLSUtil.setCellString(WS, 4,row,"r_max"); XLSUtil.setCellString(WS,
		 * 5,row,"Max. of the K function"); XLSUtil.setCellString(WS,
		 * 6,row,"p-value"); XLSUtil.setCellString(WS, 7,row,"log_10(p-value)");
		 * XLSUtil.setCellString(WS,
		 * 8,row,"Percentage of detections 2 coloc. with detections 1");
		 * XLSUtil.setCellString(WS,
		 * 9,row,"Number of detections 2 coloc. with detections 1");
		 * XLSUtil.setCellString(WS, 10,row,"Distance of Colocalization");
		 * XLSUtil.setCellString(WS, 11,row,"Coloc. Dist. Std"); row++;}
		 */
		// create the header row

		Row header = sheet.getRow(0);
		if (header == null) {
			header = sheet.createRow(0);
			header.getCell(0).setCellValue("Sequence Name");
			header.getCell(1).setCellValue("Time");
			header.getCell(2).setCellValue("Nb detections 1");
			header.getCell(3).setCellValue("Nb detections 2");
			header.getCell(4).setCellValue("r 1");
			header.getCell(5).setCellValue("r max");
			header.getCell(6).setCellValue("Max of the K function");
			header.getCell(7).setCellValue("p value");
			header.getCell(8).setCellValue("log_10(p value)");
			header.getCell(9).setCellValue("percentage of coupling");
			header.getCell(10).setCellValue("number of couples");
			header.getCell(11).setCellValue("Coupling distance");
			header.getCell(12).setCellValue("Std. of coupling distance");
		}

		// définition du ratio Z/X
		double ratio_zx = 1.0;
		if (sequence.getSizeZ() > 1) {
			ratio_zx = sequence.getPixelSizeZ() / sequence.getPixelSizeX();
		}
		// initialisation boucle en temps
		// temps total de la sequence
		T = sequence.getSizeT();
		betaCorrection(maxdist / 10, 100);
		/////////////////////////////////////////////
		///////////////////////////////////////////////////
		// initialisation
		Vector<Spot> detection_m1 = detect_.getDetectionsAtT(-1);
		Vector<Spot> detection2_m1 = detect2_.getDetectionsAtT(-1);
		Vector<Spot> detection_0 = detect_.getDetectionsAtT(0);
		Vector<Spot> detection2_0 = detect2_.getDetectionsAtT(0);
		// on construit emsuite la roi au temps t=0
		roi_t = ROI_t(dim, list_roi, 0, -1);
		// et on filtre les détections appartenant effectivement à la ROI (t)
		Vector<Spot> detection = detectionsInRoi(detection_0, detection_m1, roi_t, 0);
		Vector<Spot> detection2 = detectionsInRoi(detection2_0, detection2_m1, roi_t, 0);

		double volume = roi_t.getNumberOfPoints() * ratio_zx;
		if (roi_t.getBounds5D().isInfiniteZ()) {
			volume = volume * sequence.getSizeZ() * ratio_zx;
		}

		int nbdeta = detection.size();
		int nbdetb = detection2.size();
		// construction du distance_tab pour calculer la p_value de coloc
		double coeff_radius = 1;
		double step_min = maxdist / 10;

		if (sequence.getSizeZ() > 1) {
			mindist = Math.max(Math.pow(coeff_radius * volume / (nbdeta * nbdetb), 0.333), step_min);
		} else {
			mindist = Math.max(Math.pow(coeff_radius * volume / (nbdeta * nbdetb), 0.5), step_min);
		}
		if (mindist > maxdist) {
			mindist = maxdist;
			new AnnounceFrame("Number of spots should be insufficient for the statistical analysis");
		}
		distance_fit.clear();
		distance_fit.add((double) 0);
		distance_fit.add(mindist);
		double temp = mindist;
		if (sequence.getSizeZ() > 1) {
			while (Math.max(Math.pow(Math.pow(temp, 3) + coeff_radius * volume / (nbdeta * nbdetb), 0.333),
					temp + step_min) < maxdist) {
				temp = Math.max(Math.pow(Math.pow(temp, 3) + coeff_radius * volume / (nbdeta * nbdetb), 0.333),
						temp + step_min);
				distance_fit.add(temp);
			}
		} else {
			while (Math.max(Math.sqrt(Math.pow(temp, 2) + coeff_radius * volume / (nbdeta * nbdetb)),
					temp + step_min) < maxdist) {
				temp = Math.max(Math.sqrt(Math.pow(temp, 2) + coeff_radius * volume / (nbdeta * nbdetb)),
						temp + step_min);
				distance_fit.add(temp);
			}
		}
		int N_fit = distance_fit.size();
		////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////

		// stockage de l'histogram des distances de coloc
		double[][] histo_dist = new double[T][N_fit - 1];
		// démarrage boucle en temps
		for (int t = 0; t < T; t += 1) {
			Row row = sheet.createRow(t + 1);
			/*
			 * if (exportExcel.getValue()){ XLSUtil.setCellNumber(WS, 0, row,
			 * t);}
			 */
			// initialisation
			detection_0 = detect_.getDetectionsAtT(t);
			detection2_0 = detect2_.getDetectionsAtT(t);
			roi_t = ROI_t(dim, list_roi, t, -1);
			detection = detectionsInRoi(detection_0, detection_m1, roi_t, t);
			detection2 = detectionsInRoi(detection2_0, detection2_m1, roi_t, t);

			// calcul des parametres globaux: aire totale des rois et nb de
			// detections
			// faut-refaire une fonction calcul de volume
			volume = roi_t.getNumberOfPoints() * ratio_zx;
			if (roi_t.getBounds5D().isInfiniteZ()) {
				volume = volume * sequence.getSizeZ() * ratio_zx;
			}

			nbdeta = detection.size();
			nbdetb = detection2.size();

			////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////
			double[] delta_K = new double[N_fit - 1];
			double[][] Ktemp = new double[N_fit - 1][3];
			double[] coeffs = new double[3];
			double[] var = new double[N_fit];
			double s = -10000;
			////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////

			if (sequence.getSizeZ() > 1) {
				Ktemp = Ripley3D.correlation(sequence.getSizeZ(), roi_t, detection, detection2, volume, nbdeta, nbdetb,
						N_fit, distance_fit, ratio_zx);
			} else {
				Ktemp = Ripley2D.correlation(roi_t, detection, detection2, volume, nbdeta, nbdetb, N_fit, distance_fit);
			}

			for (int k = 0; k < N_fit - 1; k++) {
				delta_K[k] = Ktemp[k][0];
			}

			double[] retour = new double[5];
			double[] Kfit = new double[N_fit];
			Kfit[0] = 0;
			for (int i = 0; i < N_fit - 1; i++) {
				Kfit[i + 1] = Kfit[i] + delta_K[i];
			}

			if (sequence.getSizeZ() > 1) {
				var = Ripley3D.variance_theo(sequence.getSizeZ(), roi_t, detection, detection2, volume, nbdeta, nbdetb,
						N_fit, distance_fit, N_h, ratio_zx);
			} else
			// il faut calculer la roi au temps t//construire une array list
			{
				var = Ripley2D.variance_theo(roi_t, detection, detection2, volume, nbdeta, nbdetb, N_fit, distance_fit,
						N_h, results);
			}
			K.clear();
			K.add(0, 0.0);
			// calcul du sup de la fonction K normalisée

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

			double normcdf = 0.5 * (1 + ErrorFunction.erf(s / Math.sqrt(2)));
			pvalue.setValue(1 - normcdf);
			// fit sur K_fit moyen

			double[] start_est = new double[3];
			start_est[0] = 0.1D; // initial estimate of alpha
			start_est[1] = 1.01D; // initial estimate of mu
			start_est[2] = 0.3D;// sigma_ini.getValue(); // initial estimate of
								// sigma

			// voir pour entrer u
			// double sigma_fit=0.3;
			// fit_data.main(distance_fit, Kfit, area, nbdeta, start_est);
			if (sequence.getBounds5D().getSizeZ() > 1) {
				fit_data.main_manual(distance_fit, Kfit, volume, nbdeta);
				coeffs = fit_data.coeffs;
			} else {
				fit_data.main_manual(distance_fit, Kfit, volume, nbdeta);
				coeffs = fit_data.coeffs;
			}

			alpha_fit.setValue(coeffs[0]);
			mu_fit.setValue(coeffs[1]);
			number_fit.setValue(Math.ceil(alpha_fit.getValue() * nbdetb));
			sigma_fit.setValue(coeffs[2]);
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			if (export_colocalized_rois.getValue()) {

				// Une fois que l'on connait le pourcentage de coloc, on peut
				// construire une fonction qui extrait les ROIs correspondanst
				// aux spots coloc
				// 1 fonction pour déterminer la liste des apparatedSpots
				ArrayList<apparatedSpots> liste_app_detect = new ArrayList<apparatedSpots>();
				liste_app_detect = apparatedSpots.appDetectConstruction(detection, detection2);
				// 1 fonction pour construire la listes des detections 1 et 2
				// qui colocalisent à partir de la liste des apparatedDetection
				// et du pourcentage de coloc
				int ind_max = (int) (Math.ceil(alpha_fit.getValue() * nbdetb));
				ArrayList<apparatedSpots> liste_short = new ArrayList<apparatedSpots>();
				liste_short = apparatedSpots.appDetectSelect(liste_app_detect, ind_max);

				if (sequence.getSizeZ() > 1) {
					apparatedSpots.roiColoc3D(t, detections1.getValue(), detections2.getValue(), liste_short);
				} else {
					apparatedSpots.roiColoc2D(t, detections1.getValue(), detections2.getValue(), liste_short);
				}

			}

			/*
			 * if (exportExcel.getValue()) { XLSUtil.setCellNumber(WS, 1, row,
			 * nbdeta); XLSUtil.setCellNumber(WS, 2, row, nbdetb);
			 * XLSUtil.setCellNumber(WS, 3, row, mindist);
			 * XLSUtil.setCellNumber(WS, 4, row, maxdist);
			 * XLSUtil.setCellNumber(WS, 5, row, s); XLSUtil.setCellNumber(WS,
			 * 6, row, pvalue.getValue()); if (non_param.getValue()){ if
			 * (retour[4]<4){ XLSUtil.setCellNumber(WS, 7, row,
			 * Math.log10(pvalue.getValue()));} else { double x =
			 * retour[4]/Math.sqrt(2);double y =
			 * (N_fit-1)/(2*Math.sqrt(Math.PI)*x); double exponent =
			 * (Math.log(y)-Math.pow(x, 2))/Math.log(10);
			 * XLSUtil.setCellNumber(WS, 7, row, exponent); } } else{ if (s<4){
			 * XLSUtil.setCellNumber(WS, 7, row,
			 * Math.log10(pvalue.getValue()));} else { double x =
			 * s/Math.sqrt(2);double y = (1)/(2*Math.sqrt(Math.PI)*x); double
			 * exponent = (Math.log(y)-Math.pow(x, 2))/Math.log(10);
			 * XLSUtil.setCellNumber(WS, 7, row, exponent); } }
			 * XLSUtil.setCellNumber(WS, 8, row, alpha_fit.getValue());
			 * XLSUtil.setCellNumber(WS, 9, row, number_fit.getValue());
			 * XLSUtil.setCellNumber(WS, 10, row, mu_fit.getValue());}
			 * XLSUtil.setCellNumber(WS, 11, row, sigma_fit.getValue()); row++;
			 * }
			 * 
			 * /////////////////////////////// if (exportExcel.getValue()) {try
			 * { XLSUtil.saveAndClose(WW); } catch (WriteException e) { // TODO
			 * Auto-generated catch block e.printStackTrace(); } catch
			 * (IOException e) { // TODO Auto-generated catch block
			 * e.printStackTrace(); }}
			 */
			String dataSetName = "unknown sequence";
			if (input_sequence1.getValue() != null) {
				dataSetName = getDataSetName(input_sequence1.getValue());
			}
			row.getCell(0).setCellValue(dataSetName);
			row.getCell(1).setCellValue(t);
			row.getCell(2).setCellValue(nbdeta);
			row.getCell(3).setCellValue(nbdetb);
			row.getCell(4).setCellValue(mindist);
			row.getCell(5).setCellValue(maxdist);
			row.getCell(6).setCellValue(s);
			row.getCell(7).setCellValue(pvalue.getValue());
			if (s < 4) {
				row.getCell(8).setCellValue(Math.log10(pvalue.getValue()));
			} else {
				double x = s / Math.sqrt(2);
				double y = (1) / (2 * Math.sqrt(Math.PI) * x);
				double exponent = (Math.log(y) - Math.pow(x, 2)) / Math.log(10);
				row.getCell(8).setCellValue(exponent);
			}
			row.getCell(9).setCellValue(alpha_fit.getValue());
			row.getCell(10).setCellValue(number_fit.getValue());
			row.getCell(11).setCellValue(mu_fit.getValue());
			row.getCell(12).setCellValue(sigma_fit.getValue());
		}

		// export des roi qui colocalisent
		ROI[] c1 = new ROI[apparatedSpots.coloc1.size()];
		ROI[] c2 = new ROI[apparatedSpots.coloc2.size()];
		int i = 0;
		for (ROI r : apparatedSpots.coloc1) {
			c1[i] = r;
			i++;
		}
		i = 0;
		for (ROI r : apparatedSpots.coloc2) {
			c2[i] = r;
			i++;
		}

		detections1Coloc.setValue(c1);
		detections1Coloc.setValue(c2);
	}

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

	private Vector<Spot> detectionsInRoi(Vector<Spot> detection, Vector<Spot> detection_m1, ROI roi, int t) {
		// SEQUENCE A remplacer par w,h

		// create a hashMap with the detections binded to ROI

		Vector<Spot> ROIDetection = new Vector<Spot>();
		ArrayList<Spot> TestList = new ArrayList<Spot>();
		ArrayList<Spot> TestList_m1 = new ArrayList<Spot>();

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
		// fill hashMap
		for (Spot spot : detection_m1) {
			boolean b = false;
			for (Spot spot2 : TestList_m1) {
				if ((spot2.mass_center.x == spot.mass_center.x) && (spot2.mass_center.y == spot.mass_center.y)
						&& (spot2.mass_center.z == spot.mass_center.z)) {
					b = true;
				}
			}
			if (b == false)// si le spot n'a pas encore été pris en compte
			{
				TestList_m1.add(spot);
				// if (roi.contains(spot.mass_center.x, spot.mass_center.y,
				// spot.mass_center.z, t, roi.getPosition5D().getC()))
				if (roi.contains(spot.mass_center.x, spot.mass_center.y, spot.mass_center.z, roi.getPosition5D().getT(),
						roi.getPosition5D().getC())) {
					ROIDetection.add(spot);
				}
			}
		}

		if ((detection.size() == 0) && (detection_m1.size() == 0)) {
			new AnnounceFrame("There is no detection associated with the ROI(s)");

		}

		return (ROIDetection);
	}

	private ROI ROI_t(int dim, ArrayList<ROI> list_roi, int t, int c) {

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

	@Override
	public String getMainPluginClassName() {
		// TODO Auto-generated method stub
		return ColocalizationStudio.class.getName();
	}

	private String getDataSetName(Sequence sequence) {
		String dataSetName = "";
		// replace the sheet name by the file or sequence name
		dataSetName = FileUtil.getFileName(sequence.getFilename());
		if (dataSetName.isEmpty())
			dataSetName = sequence.getName();

		// make the name "safe"
		return WorkbookUtil.createSafeSheetName(dataSetName);
	}

	/**
	 * A public function which sets key variables before calling run(). Allows
	 * for easy use of plugin within a script.
	 *
	 * @param input_sequence1_inp
	 *            first input data sequence
	 * 
	 * @param input_sequence1_inp
	 *            second input data sequence
	 * 
	 * @param detections1_inp
	 *            ROIs for object detections in the first sequence
	 * 
	 * @param detections2_inp
	 *            ROIs for object detections in the second sequence
	 * 
	 * @param ROIs_inp
	 *            analysis ROIs
	 * 
	 * @param maxdistance_in_inp
	 *            Max. Radius (in pixels)
	 * 
	 */

	public void calculateColocStats(Sequence input_sequence1_inp, Sequence input_sequence2_inp,
			ArrayList<ROI> detections1_inp, ArrayList<ROI> detections2_inp, ArrayList<ROI> ROIs_inp,
			double maxdistance_in_inp) {

		// populate variables with user defined inputs
		this.input_sequence1.setValue(input_sequence1_inp);
		this.input_sequence2.setValue(input_sequence2_inp);
		this.maxdistance_in.setValue(Double.valueOf(maxdistance_in_inp));
		this.detections1.setValue((ROI[]) detections1_inp.toArray(new ROI[detections1_inp.size()]));
		this.detections2.setValue((ROI[]) detections2_inp.toArray(new ROI[detections2_inp.size()]));
		this.ROIs.setValue((ROI[]) ROIs_inp.toArray(new ROI[this.ROIs.size()]));
		
		run();
	}

	public VarROIArray getDetections1Coloc() {
		return detections1Coloc;
	}

	public VarROIArray getDetections2Coloc() {
		return detections2Coloc;
	}

	public VarDouble getNumber_fit() {
		return number_fit;
	}

	public VarDouble getAlpha_fit() {
		return alpha_fit;
	}

	public VarDouble getMu_fit() {
		return mu_fit;
	}

	public VarDouble getPvalue() {
		return pvalue;
	}

}