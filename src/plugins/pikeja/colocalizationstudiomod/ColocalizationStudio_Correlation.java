package plugins.pikeja.colocalizationstudiomod;

import java.awt.Point;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.Vector;

import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.ss.util.WorkbookUtil;

import icy.file.FileUtil;
import icy.gui.frame.progress.AnnounceFrame;
import icy.plugin.abstract_.Plugin;
import icy.plugin.interface_.PluginBundled;
import icy.roi.BooleanMask2D;
import icy.roi.BooleanMask3D;
import icy.roi.ROI;
import icy.roi.ROI2D;
import icy.sequence.Sequence;
import icy.type.point.Point5D;
import icy.util.StringUtil;
import icy.util.XLSUtil;
import jxl.write.WritableSheet;
import jxl.write.WritableWorkbook;
import jxl.write.WriteException;
import plugins.adufour.blocks.lang.Block;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.ezplug.EzVarBoolean;
import plugins.adufour.ezplug.EzVarDouble;
import plugins.adufour.ezplug.EzVarText;
import plugins.adufour.vars.lang.VarBoolean;
import plugins.adufour.vars.lang.VarDouble;

import plugins.adufour.vars.lang.VarDoubleArrayNative;
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

public class ColocalizationStudio_Correlation extends Plugin implements Block, PluginBundled {	
				
	
	
	VarROIArray	ROIs		= new VarROIArray("Analysis ROIs");
	
	VarSequence input_sequence1 = new VarSequence("Input sequence 1", null);
	VarSequence input_sequence2 = new VarSequence("Input sequence 2", null);
	
	VarROIArray	ROIs_M1		= new VarROIArray("ROIs for Manders analysis Seq1");
	VarROIArray	ROIs_M2		= new VarROIArray("ROIs for Manders analysis Seq2");
	VarBoolean surface = new VarBoolean("Consider pixel surface for Manders analysis", false);
	VarDouble pearson_coeff = new VarDouble("Pearson coeff. R", (double) 0.0);
	VarDouble pvalue_pearson = new VarDouble("p value Pearson (Randomization)",0.0);
	
	VarDouble M1 = new VarDouble("Manders M1", 0.0);
	VarDouble pvalue_M1 = new VarDouble("p value M1 (Randomization)", 0.0);
	VarDouble M2 = new VarDouble("Manders M2", 0.0);
	VarDouble pvalue_M2 = new VarDouble("p value M2 (Randomization)", 0.0);
	
	
	VarDouble ICCS1=new VarDouble("Cross-Correlation 1", 0.0);
	VarDouble ICCS2=new VarDouble("Cross-Correlation 2", 0.0);
			
	
	/*VarBoolean                          exportExcel        = new VarBoolean("Export results in excel",false);	
	VarFile                             exportExcelFile        = new VarFile("Excel File", null);*/
	VarWorkbook                book              = new VarWorkbook("Workbook", (Workbook) null);
	
	ArrayList<ROI> list_roi = new ArrayList<ROI>();	
	@Override
	public void declareInput(VarList inputMap)
	{
		
		inputMap.add("Input Sequence 1", input_sequence1);
		inputMap.add("Input Sequence 2", input_sequence2);
		inputMap.add("Analysis ROIs", ROIs);
		inputMap.add("ROIs for Manders Seq. 1", ROIs_M1);
		inputMap.add("ROIs for Manders Seq. 2", ROIs_M2);
		inputMap.add("Consider pixel surface for Manders Analysis", surface);
		
		/*inputMap.add( "Excel Export", exportExcel );
		inputMap.add( "Excel File", exportExcelFile );*/
		}

	@Override
	public void declareOutput(VarList outputMap)
	{		
		outputMap.add("Pearson Coeff. R",pearson_coeff);
		outputMap.add("Pearson p-value",pvalue_pearson);
		outputMap.add("Manders M1",M1);
		outputMap.add("M1 p-value",pvalue_M1);
		outputMap.add("Manders M2",M2);
		outputMap.add("M2 p-value",pvalue_M2);
		outputMap.add("Cross-Correlation 1",ICCS1);
		outputMap.add("Cross-Correlation 2",ICCS2);
		outputMap.add("Workbook", book);
	}
	
	@Override
	public void run() {
		if (input_sequence1.getValue() == null || input_sequence2.getValue() == null) {
			new AnnounceFrame("Please first select sequences");
			return;} 
    	if (input_sequence1.getValue().getSizeX()!= input_sequence2.getValue().getSizeX() ||input_sequence1.getValue().getSizeY()!= input_sequence2.getValue().getSizeY() || input_sequence1.getValue().getSizeZ()!= input_sequence2.getValue().getSizeZ()||input_sequence1.getValue().getSizeT()!= input_sequence2.getValue().getSizeT()) {
			new AnnounceFrame("Sequences must have the same (X,Y,Z,T) dimensions");
			return;}
    	//créaton du workbook
		 if (book.getValue()==null) 
			 {book.setValue(new HSSFWorkbook());}
		 book.getValue().setMissingCellPolicy(Row.CREATE_NULL_AS_BLANK);  
    	performAnalysis(input_sequence1, input_sequence2);  
	}
	private void performAnalysis(VarSequence sequence1,VarSequence sequence2) {
		// initialisation xls
				/*int row = 0;		
				WritableWorkbook WW = null;
				WritableSheet WS = null;
				if (exportExcel.getValue()) {
					int page = 1;
					try {
						File f = exportExcelFile.getValue(true);
						if (!FileUtil.getFileExtension(f.getPath(), false)
								.equalsIgnoreCase("xls"))
							f = new File(f.getPath() + ".xls");
						WW = XLSUtil.loadWorkbookForWrite(f); 
						}							
					 catch (Exception e) {
						e.printStackTrace();
						return;}			
					WS  = XLSUtil.createNewPage(WW, "Page" + page);
					XLSUtil.setCellString(WS, 0, 0, "Date of XLS page:");
					row++;
					XLSUtil.setCellString(WS, 0, row, new Date().toString());			
					row++;}*/
		//initialisation du workbook
				Workbook wb = book.getValue();
		       // create the sheet
		        String sheetName = "Coloc. Object Analysis";        
		        Sheet sheet = wb.getSheet(sheetName);
		        if (sheet == null) sheet = wb.createSheet(sheetName);
		        
		        
				Sequence sequence = input_sequence1.getValue();
				int dim =2;
				if (sequence.getSizeZ()>1)
				dim = 3;
				
				//gestion des rois d'analyse en entrée		
				if (ROIs.getValue()==null||ROIs.getValue().length==0)
				{
					for (int t=0;t<input_sequence1.getValue().getSizeT();t++)
					{
						ROI roi=null;				
					
					ROI2DRectangle r= new ROI2DRectangle(sequence.getBounds2D());
					for (int h=0;h<sequence.getSizeZ();h++)
					{
						r.setZ(h);r.setT(t);						
						roi = r.getUnion(roi);
					}
					list_roi.add(roi);
					}
				}
				else
				{ 
					for (ROI r:ROIs.getValue())
					{list_roi.add(r);}
				}
					
					//on teste si les dimension en c sonts compatible
					double c1 = list_roi.get(0).getPosition5D().getC();
					for (ROI r:list_roi)
					{double c = r.getPosition5D().getC();
					if ((c!=c1) && (c!=(-1))){new AnnounceFrame("ROI channels are incompatibles");return;}
					}

					//on teste si les dimensions en temps/z sont incompatibles pour une union
					boolean one_z=false;
					boolean all_z = true;
					
					for (ROI r:list_roi)
					{	
						if (r.getBounds5D().isInfiniteZ()==false)
						{
							all_z=false;					
						}
						if (r.getBounds5D().isInfiniteZ())
						{
							one_z=true;					
						}
					}
					//gestion de l'exception
					if (one_z==true && all_z==false)
					{new AnnounceFrame("Incompatibility in Z dimensions between ROIs");
					return;}	
				
				
				//on construit emsuite la roi d'analyse au temps t=0
				ROI roi_t = ROI_t(dim,list_roi,0,-1);		
				int T;
				/*if (exportExcel.getValue()) {
		    		
		    		XLSUtil.setCellString(WS, 0, row, "Time");
		    		XLSUtil.setCellString(WS, 1, row, "Pearson coefficient");
		    		XLSUtil.setCellString(WS, 2, row, "p value");
		    		XLSUtil.setCellString(WS, 3, row, "M1");
		    		XLSUtil.setCellString(WS, 4, row, "p value M1");
		    		XLSUtil.setCellString(WS, 5, row, "M2");
		    		XLSUtil.setCellString(WS, 6, row, "p value M2");
		    		XLSUtil.setCellString(WS, 7, row, "ICCS1");
		    		XLSUtil.setCellString(WS, 8, row, "ICCS2");					row++;}*/
// create the header row
		        
		        Row header = sheet.getRow(0);
		        if (header == null)
		        {
		            header = sheet.createRow(0);
		            header.getCell(0).setCellValue("Sequence Name");
		            header.getCell(1).setCellValue("Time");
		            header.getCell(2).setCellValue("Pearson coefficient");
		            header.getCell(3).setCellValue("p value");
		            header.getCell(4).setCellValue("M1");
		            header.getCell(5).setCellValue("p value (M1)");
		            header.getCell(6).setCellValue("M2");
		            header.getCell(7).setCellValue("p value (M2)");
		            header.getCell(8).setCellValue("Cross Corr. 1");
		            header.getCell(9).setCellValue("Cross Corr. 2");		            
		        }	
											
			//initialisation boucle en temps
				// temps total de la sequence		
				T=input_sequence1.getValue().getSizeT();
					//démarrage boucle en temps
					for (int t = 0; t < T; t += 1) {
						Row row = sheet.createRow(t+1);
						/*if (exportExcel.getValue()){
							XLSUtil.setCellNumber(WS, 0,row,t);	}*/						
						
						roi_t = ROI_t(dim,list_roi,t,-1);
						
					double pearson[] = Correlation.pearson_TCL(input_sequence1.getValue(),input_sequence2.getValue(),t,0,0,roi_t);
					pearson_coeff.setValue(pearson[0]);					
					pvalue_pearson.setValue(pearson[1]);
					
					/*//gerer l'export excel
					if (exportExcel.getValue()) {
						XLSUtil.setCellNumber(WS, 1,row,pearson_coeff.getValue()); 
						XLSUtil.setCellNumber(WS, 2,row,pvalue_pearson.getValue());}*/
					String dataSetName = "unknown sequence";		 
		            if (input_sequence1.getValue()!=null){dataSetName = getDataSetName(input_sequence1.getValue());} 
		            row.getCell(0).setCellValue(dataSetName);    	
		            row.getCell(1).setCellValue(t);			
		            row.getCell(2).setCellValue(pearson_coeff.getValue());
		            row.getCell(3).setCellValue(pvalue_pearson.getValue());
					
					//on enleve ensuite les ROIs "d'étude"
					ArrayList<ROI> list1_=new ArrayList<ROI>();
					ArrayList<ROI> list2_=new ArrayList<ROI>();
					for (ROI r:ROIs_M1)
					{if (list_roi.contains(r)){}else{list1_.add(r);}}
					for (ROI r:ROIs_M2)
					{if (list_roi.contains(r)){}else{list2_.add(r);}}
					//on trie les ROIs avec le bon canal
					//on trie les ROIs avec le bon canal
					ArrayList<ROI> list1__=new ArrayList<ROI>();
					ArrayList<ROI> list2__=new ArrayList<ROI>();
					for (ROI r:list1_)
					{if ((r.getPosition5D().getC()==0)||(r.getBounds5D().isInfiniteC())){list1__.add(r);}}
					for (ROI r:list2_)
					{if ((r.getPosition5D().getC()==0)||(r.getBounds5D().isInfiniteC())){list2__.add(r);}}
					//
					
					ROI roi1 = ROI_t(dim,list1__,t,0);
					ROI roi2 = ROI_t(dim,list2__,t,0);
					double manders[] = new double[4];
					if (roi1==null || roi2==null)
					{}
					else{
					manders = Correlation.MandersCoeff(input_sequence1.getValue(), input_sequence2.getValue(), t, 0, 0,surface.getValue(), roi1, roi2,roi_t);
					M1.setValue(manders[0]);
					M2.setValue(manders[1]);
								
					pvalue_M1.setValue(manders[2]);
					pvalue_M2.setValue(manders[3]);}
					
					/*if (exportExcel.getValue()) {
						XLSUtil.setCellNumber(WS, 3,row,M1.getValue()); 
						XLSUtil.setCellNumber(WS, 4,row,pvalue_M1.getValue()); 
						XLSUtil.setCellNumber(WS, 5,row,M2.getValue()); 
						XLSUtil.setCellNumber(WS, 6,row,pvalue_M2.getValue()); }}*/
					row.getCell(4).setCellValue(M1.getValue());
		            row.getCell(5).setCellValue(pvalue_M1.getValue());
		            row.getCell(6).setCellValue(M2.getValue());
		            row.getCell(7).setCellValue(pvalue_M2.getValue());
		            
					double results_iccs[] = Correlation.ICCS_compute(input_sequence1.getValue(),input_sequence2.getValue(),0,0,t,roi_t);
						ICCS1.setValue(results_iccs[0]);
						ICCS2.setValue(results_iccs[1]);
						/*if (exportExcel.getValue()) {
							XLSUtil.setCellNumber(WS, 7,row,ICCS1.getValue()); 
							XLSUtil.setCellNumber(WS, 8,row,ICCS2.getValue()); row++;}*/
						row.getCell(8).setCellValue(ICCS1.getValue());
			            row.getCell(9).setCellValue(ICCS2.getValue());			            
					
					}					
	}
								
				/*if (exportExcel.getValue()) {
					try {
						XLSUtil.saveAndClose(WW);
					} catch (WriteException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}}*/

		 		

	
			
		
			public Point2D getMassCenter(ROI2D roi)
		   {
		       double x = 0, y = 0;
		       long len = 0;			       

		       final BooleanMask2D mask = roi.getBooleanMask(true);
		       final boolean m[] = mask.mask;
		       final int h = mask.bounds.height;
		       final int w = mask.bounds.width;

		       int off = 0;
		       for (int j = 0; j < h; j++)
		       {
		           for (int i = 0; i < w; i++)
		           {
		               if (m[off++])
		               {
		                   x += i;
		                   y += j;
		                   len++;
		               }
		           }
		       }

		       final Point2D pos2d = roi.getPosition2D();
		       return new Point2D.Double(pos2d.getX() + (x / len), pos2d.getY() + (y / len));
		   }

	
			private ROI ROI_t(int dim,ArrayList<ROI> list_roi,int t,int c)
			{
				
				//check compatibility en z si dim =3
						if (dim==3){
						boolean one_z=false;
						boolean all_z = true;		
						for (ROI r0:list_roi)
						{	
							if (r0.getBounds5D().isInfiniteZ()==false)
							{
								all_z=false;					
							}
							if (r0.getBounds5D().isInfiniteZ())
							{
								one_z=true;					
							}
						}
						//gestion de l'exception
						if (one_z==true && all_z==false)
						{return null;}}	
						
				ROI roi_t=null;
				
				for (ROI r:list_roi)
				{
					Point5D pt0 = r.getPosition5D();
					if (dim==2){
						pt0.setZ(-1);}
					if ((r.getBounds5D().isInfiniteC())||(c==-1)){pt0.setC(c);}		
					if ((r.getBounds5D().isInfiniteT())||(t==-1)){pt0.setT(t);}
					r.setPosition5D(pt0);
					if ((pt0.getT()==t) && (pt0.getC()==c))
					{
						roi_t=r.getUnion(roi_t);
					}
				}
				return roi_t;
			}

	@Override
	public String getMainPluginClassName() {
		// TODO Auto-generated method stub
		return ColocalizationStudio.class.getName();
	}
	private String getDataSetName(Sequence sequence)
    {
        String dataSetName = "";
            // replace the sheet name by the file or sequence name
            dataSetName = FileUtil.getFileName(sequence.getFilename());
            if (dataSetName.isEmpty()) dataSetName = sequence.getName();
        
        // make the name "safe"
        return WorkbookUtil.createSafeSheetName(dataSetName);
    }
	
}