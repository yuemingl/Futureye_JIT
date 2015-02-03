package edu.uta.futureye.application;

import edu.uta.futureye.algebra.SparseVectorHashMap;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.application.HumanReal.Part;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.function.operator.FMath;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.util.DoubleRange;
import edu.uta.futureye.util.container.NodeList;

public class HumanRealPost {
	/**
	 * @param args
	 */	
	public static void main(String[] args) {
//		computeHboHb(Part.LEFT,"mesh3DLeft.grd",
//				"138-average/750_tail2",
//				"138-average/850_tail2",
//				"138-average/Hb_tail2");
//		computeHboHb(Part.RIGHT,"mesh3DRight.grd",
//				"138-average/750_tail2",
//				"138-average/850_tail2",
//				"138-average/Hb_tail2");
		computeHboHb(Part.RIGHT,"mesh3DRight.grd",
				"138-average/750avg",
				"138-average/850avg",
				"138-average/Hbavg");		
		computeHboHb(Part.LEFT,"mesh3DLeft.grd",
						"138-average/750avg",
						"138-average/850avg",
						"138-average/Hbavg");

	}

	/**
	 * 计算Hbo(含氧血红素)和Hb(缺氧血红素)
	 * 
	 * (OD1) = (EHbo1 EHb1) * (Hbo)
	 * (OD2)   (EHbo2 EHb2)   (Hb )
	 * 
	 * =>
	 * 
	 * (Hbo) = inv[(EHbo1 EHb1)] * (OD1)
	 * (Hb )      [(EHbo2 EHb2)]   (OD2)
	 * 
	 * OD1: 750nm 重构的吸收系数mu_a (注：OD(Optical Density) = absorbance = mu_a != 3*mu_a*mu_s')
	 * OD2: 850nm 重构的吸收系数mu_a
	 * 
	 */
	public static void computeHboHb(Part part,String gridName, 
			String wave750Folder, String wave850Folder,
			String outFolder) {
		String baseFolder = "./HumanReal/Output/";
		//3D Mesh
		MeshReader reader = new MeshReader("./HumanReal/"+gridName);
		Mesh mesh = reader.read3DMesh();
		mesh.computeNodeBelongsToElements();
		mesh.computeNeighborNodes();
		//GetExtinctions([750 850])
		//Mat=[ 1381.800   3528.196;
		//      2526.391   1798.643 ]
		double[] invMat ={
				 -0.000279803435628786,   0.000548858980004226,
				  0.000393014556830702,  -0.000214957825066929
		   };
		
		String format750 = null;
		String format850 = null;
		if(part==Part.LEFT) {
			format750 = baseFolder+wave750Folder+"/3D_t%05d_Left.dat";
			format850 = baseFolder+wave850Folder+"/3D_t%05d_Left.dat";
		} else {
			format750 = baseFolder+wave750Folder+"/3D_t%05d_Right.dat";
			format850 = baseFolder+wave850Folder+"/3D_t%05d_Right.dat";
		}
		String file750,file850;
		Vector v750 = null;
		Vector v850 = null;
		
		int N=25;
		Vector Hb=null,Hbo=null,HbT=null;
		Vector HbPeak=new SparseVectorHashMap(N);
		Vector HboPeak=new SparseVectorHashMap(N);
		Vector HbTPeak=new SparseVectorHashMap(N);
		Vector index=new SparseVectorHashMap(N);
		DoubleRange rx = new DoubleRange(4.0,6.5);
		DoubleRange ry = new DoubleRange(-5.5,-2.5);
		DoubleRange rz = null;
		if(part==Part.LEFT) {
			rz = new DoubleRange(3.5,4.5);
		} else {
			rz = new DoubleRange(12.2,13.8);
		}
		for(int i=0;i<N;i++) {
			file750 = String.format(format750,i);
			file850 = String.format(format850,i);
			
			v750 = DataReader.readVector(file750,4).scale(300);
			v850 = DataReader.readVector(file850,4).scale(300);
			cutByRange(mesh,v750,rx,ry,rz);
			cutByRange(mesh,v850,rx,ry,rz);
			Tools.plotVector(mesh, baseFolder+outFolder, String.format("750cut_%05d"+part+".dat",i), v750);
			Tools.plotVector(mesh, baseFolder+outFolder, String.format("850cut_%05d"+part+".dat",i), v850);
			
			Hbo = v750.copy();
			Hbo.scale(invMat[0]).add(invMat[1], v850);
			Hb = v750.copy();
			Hb.scale(invMat[2]).add(invMat[3], v850);
			
//			Hb = Utils.gaussSmooth(mesh, Hb, 2, 0.5);
//			Hb = Utils.gaussSmooth(mesh, Hb, 2, 0.5);
//			Hbo = Utils.gaussSmooth(mesh, Hbo, 2, 0.5);
//			Hbo = Utils.gaussSmooth(mesh, Hbo, 2, 0.5);
			
			HbT = FMath.axpy(1.0, Hb, Hbo);
			
			Node maxNode = new Node();
			//固定点Hbo的最大点
			if(part==Part.LEFT) {
				HbPeak.set(i+1,localMax(mesh,Hb,rz,maxNode));
			} else {
				HbPeak.set(i+1,localMax(mesh,Hb,rz,maxNode));
			}
			HboPeak.set(i+1,Hbo.get(maxNode.getIndex()));
			HbTPeak.set(i+1,HbT.get(maxNode.getIndex()));
			index.set(i+1,i-4);
			
			Tools.plotVector(mesh, baseFolder+outFolder, String.format("Hb_%05d"+part+".dat",i), Hb);
			Tools.plotVector(mesh, baseFolder+outFolder, String.format("Hbo_%05d"+part+".dat",i), Hbo);
			Tools.plotVector(mesh, baseFolder+outFolder, String.format("HbT_%05d"+part+".dat",i), HbT);
			
		}
		
		DataWriter.writeVector(baseFolder+outFolder+"/HbPeak"+part+".dat", index, HbPeak);
		DataWriter.writeVector(baseFolder+outFolder+"/HboPeak"+part+".dat", index, HboPeak);
		DataWriter.writeVector(baseFolder+outFolder+"/HbTPeak"+part+".dat", index, HbTPeak);
	}
	
	public static double localMax(Mesh mesh, Vector v, DoubleRange rz, Node maxNode) {
		double max = -1e17;
		NodeList nodes = mesh.getNodeList();
		int maxIndex = 0;
		for(int i=1;i<v.getDim();i++) {
			double z = nodes.at(i).coord(3);
			if(rz.isIncLR(z)) {
				double val = v.get(i);
				if(val > max) {
					maxIndex = i;
					max  = val;
				}
			}
		}
		if(maxNode != null) {
			maxNode.set(maxIndex, nodes.at(maxIndex).coords());
		}
		return max;
	}
	
	public static void cutByRange(Mesh mesh, Vector v, DoubleRange rx, DoubleRange ry, DoubleRange rz) {
		NodeList nodes = mesh.getNodeList();
		for(int i=1;i<v.getDim();i++) {
			double x = nodes.at(i).coord(1);
			double y = nodes.at(i).coord(2);
			double z = nodes.at(i).coord(3);
			if(rx.isIncLR(x) && ry.isIncLR(y) && rz.isIncLR(z)) {
				//System.out.println("In  x="+x+" y="+y+" z="+z);
			} else {
				v.set(i, 0.0);
				//System.out.println("Out x="+x+" y="+y+" z="+z);
			}
		}
	}

}
