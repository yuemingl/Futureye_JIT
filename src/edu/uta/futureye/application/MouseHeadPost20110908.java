package edu.uta.futureye.application;

import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.function.FMath;
import edu.uta.futureye.io.MeshReader;

public class MouseHeadPost20110908 {
	/**
	 * @param args
	 */	
	public static void main(String[] args) {
		computeHboHb("mouse0908_1");
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

		String format = ".\\MouseHead20110908\\Results\\830nm_"+gridName+
						"\\%s\\alpha_omega_BL%03d_ext_smooth.dat";
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
//			Tools.plotVector(meshOmega, "MouseHead20110908\\Results\\760nm\\Average", 
//					String.format("avg%03d.dat",i), v1);
			Tools.plotVector(meshOmega, "MouseHead20110908\\Results\\830nm\\Average", 
					String.format("avg%03d.dat",i), v1);
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
	 * OD1: 760nm 重构的吸收系数mu_a (注：OD(Optical Density) = absorbance = mu_a != 3*mu_a*mu_s')
	 * OD2: 830nm 重构的吸收系数mu_a
	 * 
	 * @param gridName "mouse","mouse1",...,"mouse5"
	 */
	public static void computeHboHb(String gridName) {
		//老鼠头，梯形区域
		MeshReader reader = new MeshReader(gridName+"_omega.grd");
		Mesh meshOmega = reader.read2DMesh();
		//Mat = [1486.5865 3843.707; 2321.424 1791.734];
		double[] invMat ={
			   -0.000286251218678230,  0.000614078771174764,
				0.000370875614945573, -0.000237500207785087
		   };
		
		//手工建立一下文件夹，并复制结果文件
		//String folderPostfix = "_factor=3000_rat2_S13S2";
		//String folderPostfix = "_baseline_rat1_S9S2";
		//String folderPostfix = "_baseline_rat1_S17S2ok";
		//String folderPostfix = "_baseline_rat1_S9S2ok";
		
		//String folderPostfix = "_baseline_rat1_S13S2ok";
		//String folderPostfix = "_baseline_rat2_S13S2ok";
		
		//String folderPostfix = "_baseline_rat2_S13S2_99ok";
		//String folderPostfix = "_baseline_rat2_S13S2_99mua=0.1";
		String folderPostfix = "_baseline_rat2_S13S2_341mua=0.1";
		//String folderPostfix = "_baseline_rat1_S13S2_261mua=0.1";
		
		String format760 = ".\\MouseHead20110908\\Results\\760nm_"+gridName+folderPostfix+
						"\\final_alpha_omega_BL%03d.dat";
		String format830 = ".\\MouseHead20110908\\Results\\830nm_"+gridName+folderPostfix+
						"\\final_alpha_omega_BL%03d.dat";
		String file760,file830;
		Vector v1 = null;
		Vector v2 = null;
		
		
		//int [] set = {20,30};
		//for(int i:set) {
		//for(int i=1;i<=37;i++) {
		//for(int i=1;i<=31;i++) {
		int N=341;
		
		Vector Hb=null,Hbo=null,HbT=null;
		for(int i=2;i<=N;i++) {
			file760 = String.format(format760,i);
			file830 = String.format(format830,i);
//rat1			
//			v1 = DataReader.readVector(file760).shift(0.085);//.scale(10);//shift
//			v2 = DataReader.readVector(file830).shift(0.085);//.scale(10);//shift
//rat2			
			v1 = DataReader.readVector(file760).shift(0.2);//.scale(10);//shift
			v2 = DataReader.readVector(file830).shift(0.2);//.scale(10);//shift
			Hbo = v1.copy();
			Hbo.scale(invMat[0]).add(invMat[1], v2);
			Hb = v1.copy();
			Hb.scale(invMat[2]).add(invMat[3], v2);
			HbT = FMath.axpy(1.0, Hb, Hbo);

			Tools.plotVector(meshOmega, "MouseHead20110908\\Results\\"+gridName+folderPostfix,
					String.format("Hbo%03d.dat",i), Hbo);
			Tools.plotVector(meshOmega, "MouseHead20110908\\Results\\"+gridName+folderPostfix,
					String.format("Hb%03d.dat",i), Hb);
			Tools.plotVector(meshOmega, "MouseHead20110908\\Results\\"+gridName+folderPostfix,
					String.format("HbT%03d.dat",i), HbT);
			
//			for(int k=1;k<=Hbo.getDim();k++) {
//				if(Math.abs(Hbo.get(k)+Hb.get(k))<Constant.eps)
//					Hbo.set(k,0);
//				else					
//					Hbo.set(k,Hbo.get(k)/(Hbo.get(k)+Hb.get(k)));
//			}
			double max = FMath.max(Hbo);
			max += FMath.max(Hb); //HbT=Hbo+Hb
			Vector HboPercent =  Hbo.copy().scale(1/max);
			Tools.plotVector(meshOmega, "MouseHead20110908\\Results\\"+gridName+folderPostfix,
					String.format("HboPercent%03d.dat",i),HboPercent);
			
			//将标准输出的内容复制到_MaxHbo.dat（需要新建）文件
			//Index Hbo/(HbT)
			System.out.println(i+" "+100.0*HboPercent.get(256)+" "+100.0*HboPercent.get(249));
			
			
			//Index Hbo(R) Hb(R) Hbo(L) Hb(L)
//			int[] nNodeOnRightBrain = {186,207,228,249};//165;//251;
//			int[] nNodeOnLeftBrain = {193,214,235,256};//299;//256;
//			double HboR=0.0,HbR=0.0,HboL=0.0,HbL=0.0,HbTR=0.0,HbTL=0.0;
//			int avgN = nNodeOnRightBrain.length;
//			for(int k=0;k<avgN;k++) {
//				HboR += Hbo.get(nNodeOnRightBrain[k]);
//				HbR  += Hb.get(nNodeOnRightBrain[k]);
//				HboL += Hbo.get(nNodeOnLeftBrain[k]);
//				HbL  += Hb.get(nNodeOnLeftBrain[k]);
//			}
//			HboR /= avgN;
//			HbR /= avgN;
//			HboL /= avgN;
//			HbL /= avgN;
//			HbTR = HboR+HbR;
//			HbTL = HboL+HbL;
//			
//			System.out.println(String.format("%d %.8f %.8f %.8f %.8f %.8f %.8f", 
//				i,HboR,HbR,HboL,HbL,HbTR,HbTL));			
		}
	}

}
