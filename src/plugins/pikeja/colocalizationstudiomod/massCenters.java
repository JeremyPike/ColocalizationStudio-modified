package plugins.pikeja.colocalizationstudiomod;

import icy.roi.BooleanMask2D;
import icy.roi.ROI;
import icy.roi.ROI2D;
import icy.type.point.Point3D;

import java.awt.geom.Point2D;

import plugins.kernel.roi.roi2d.ROI2DPoint;

public class massCenters {
	public static Point2D getMassCenter2D(ROI2D roi)
	   {
		if (roi instanceof ROI2DPoint)
		{
			return new Point2D.Double(roi.getBounds2D().getX(),roi.getBounds2D().getY());
		}
		else
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
	   }
	public static Point3D getMassCenter(ROI roi)
	   {
		if (roi instanceof ROI2DPoint) {							
			return new Point3D.Double(roi.getPosition5D().getX(),roi.getPosition5D().getY(),roi.getPosition5D().getZ());				
		}
		
	       double x = 0, y = 0, z = 0;
	       long len = 0;			       
	       final double t = roi.getPosition5D().getT();
	       final double c = roi.getPosition5D().getC();
	       final int minX =(int) roi.getBounds5D().getMinX();
	       final int maxX = (int)roi.getBounds5D().getMaxX();
	       final int minY =(int) roi.getBounds5D().getMinY();
	       final int maxY = (int)roi.getBounds5D().getMaxY();
	       final int minZ, maxZ;
	       
	       if (roi.getBounds5D().isInfiniteZ())
	       {for (int j = minY; j < maxY+1; j++)
	       {
	           for (int i = minX; i < maxX+1; i++)
	           {		        	   
	               if (roi.contains(i, j, -1,t , c))
	               {
	                   x += i;
	                   y += j;		                   
	                   len++;
	               }
	           }
	       }		       
	       x/=len;y/=len;z=0;
	       }
	       else
	       {
	       minZ =(int) roi.getBounds5D().getMinZ();
	       maxZ = (int)roi.getBounds5D().getMaxZ();		       		       
	       for (int k = minZ; k < maxZ; k++){
	       for (int j = minY; j < maxY; j++)
	       {
	           for (int i = minX; i < maxX; i++)
	           {		        
	               if (roi.contains(i, j, k, t, c))
	               {
	                   x += i;
	                   y += j;
	                   z+= k;
	                   len++;
	               }
	           }
	       }}		       
	       x/=len;y/=len;z/=len;}		       
	       return new Point3D.Double(x,y,z);
	   }

}


