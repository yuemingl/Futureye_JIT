package edu.uta.futureye.application;

import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.function.operator.FMath;
import edu.uta.futureye.io.MeshReader;

public class MouseHeadPost {
	/**
	 * @param args
	 */	
	public static void main(String[] args) {
		//averageData("830nm","mouse");
		computeHboHb("mouse4");
	}

	/**
	 * 计算吸收系数mu_a的平均结果（不同的BL tail光源组合）
	 * @param waveType "760nm","830nm"
	 * @param gridName "mouse","mouse1",...,"mouse5"
	 */
	public static void averageData(String waveType, String gridName) {
		//老鼠头，梯形区域
		MeshReader reader = new MeshReader(gridName+"_omega.grd");
		Mesh meshOmega = reader.read2DMesh();
		
		//////////////////////////////Configure///////////////
		//不同的BL tail光源组合
		String[] BxLx = {
				"B12L6",
				"B12L7",
				"B12L8",
				"B12L9"
		};
		///////////////////////////////////////////////////////

		String format = ".\\MouseHead\\Results\\830nm_"+gridName+
						"\\%s\\alpha_omega_BL%02d_ext_smooth.dat";
		String file;
		Vector v1 = null;
		Vector v2 = null;
		for(int i=1;i<=37;i++) {
			file = String.format(format, BxLx[0],i);
			v1 = DataReader.readVector(file);
			v2 = null;

			for(int j=1;j<BxLx.length;j++) {
				file = String.format(format, BxLx[j],i);
				v2 = DataReader.readVector(file);
				v1.add(v2);
			}
			v1.scale(1.0/BxLx.length);
//			Tools.plotVector(meshOmega, "MouseHead\\Results\\760nm\\Average", 
//					String.format("avg%02d.dat",i), v1);
			Tools.plotVector(meshOmega, "MouseHead\\Results\\830nm\\Average", 
					String.format("avg%02d.dat",i), v1);
		}
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
	 * OD1: 760nm 重构的吸收系数mu_a
	 * OD2: 830nm 重构的吸收系数mu_a
	 * 
	 * @param gridName "mouse","mouse1",...,"mouse5"
	 */
	public static void computeHboHb(String gridName) {
		//老鼠头，梯形区域
		MeshReader reader = new MeshReader(gridName+"_omega.grd");
		Mesh meshOmega = reader.read2DMesh();
		double[] invMat ={
			   -286.251218678230,  614.078771174764,
				370.875614945573, -237.500207785087
		   };
		
		String format760 = ".\\MouseHead\\Results\\760nm_"+gridName+
						"\\final_alpha_omega_BL%02d.dat";
		String format830 = ".\\MouseHead\\Results\\830nm_"+gridName+
						"\\final_alpha_omega_BL%02d.dat";
		String file760,file830;
		Vector v1 = null;
		Vector v2 = null;
		//int [] set = {20,30};
		//for(int i:set) {
		for(int i=1;i<=37;i++) {
			file760 = String.format(format760,i);
			file830 = String.format(format830,i);
			v1 = DataReader.readVector(file760);
			v2 = DataReader.readVector(file830);
			Vector Hbo = v1.copy();
			Hbo.scale(invMat[0]).add(invMat[1], v2);
			Vector Hb = v1.copy();
			Hb.scale(invMat[2]).add(invMat[3], v2);

			Tools.plotVector(meshOmega, "MouseHead\\Results\\"+gridName,
					String.format("Hbo%02d.dat",i), Hbo);
			Tools.plotVector(meshOmega, "MouseHead\\Results\\"+gridName,
					String.format("Hb%02d.dat",i), Hb);
			
//			for(int k=1;k<=Hbo.getDim();k++) {
//				if(Math.abs(Hbo.get(k)+Hb.get(k))<Constant.eps)
//					Hbo.set(k,0);
//				else					
//					Hbo.set(k,Hbo.get(k)/(Hbo.get(k)+Hb.get(k)));
//			}
			double max = FMath.max(Hbo);
			max += FMath.max(Hb);
			Tools.plotVector(meshOmega, "MouseHead\\Results\\"+gridName,
					String.format("HboPercent%02d.dat",i), Hbo.scale(1/max));
			System.out.println(i+" "+100*FMath.max(Hbo));
		}
	}

}
