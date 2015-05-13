package edu.uta.futureye.application;

import java.util.ArrayList;
import java.util.HashMap;

import edu.uta.futureye.algebra.FullMatrix;
import edu.uta.futureye.algebra.SpaceVector;
import edu.uta.futureye.algebra.SparseMatrixRowMajor;
import edu.uta.futureye.algebra.SparseVectorHashMap;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.FMath;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.Vector2Function;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.io.MatlabMatFileReader;
import edu.uta.futureye.io.MatlabMatFileWriter;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.Utils;
import edu.uta.futureye.util.container.NodeList;
import static edu.uta.futureye.function.FMath.*;

public class BOEPaper {
	String outputFolder = "BOEPaper";
	String gridFileBig = "OpticalPaperBig.grd";
	String gridFileSmall = "OpticalPaperSmall.grd";
	Mesh meshBig = null;
	Mesh meshSmall = null;
	ModelDOT model = new ModelDOT();
	//Parameters
	double init_mu_sp = 10.0; //mu_sp = 10.0 or 50.0/3.0
	double mu_sp_error = 0.05;
	Function mu_sp = null; //mu_sp
	Function mu_a = null; //mu_a
	boolean debug = true;
	int[] lightNums = {1, 2, 3, 4};
	double[][] lightPosition = {{-2.0,0.0}, {0.0, 2.0}, {0.0, -2.0}, {2.0, 0.0}};
	String[] lightPostfixs = {"left","top","bottom","right"};
	
	public BOEPaper() {
		MeshReader readerForward = new MeshReader(gridFileBig);
		meshBig = readerForward.read2DMesh();
		MeshReader readerGCM = new MeshReader(gridFileSmall);
		meshSmall = readerGCM.read2DMesh();
		
		//Use element library to assign degree of freedom (DOF) to element
		Tools.assignLinearShapFunction(meshBig);
		Tools.assignLinearShapFunction(meshSmall);
		meshBig.computeNodeBelongsToElements();
		meshBig.computeNeighborNodes();
		meshSmall.computeNodeBelongsToElements();
		meshSmall.computeNeighborNodes();
	}
	
	/**
	 * 1: light at the left
	 * 2: light at the top
	 * 3: light at the bottom
	 * 4: light at the right
	 * @param lightNum
	 */
	private void setLightPosition(int lightNum) {
		switch(lightNum) {
		case 1:
			model.setLightPosition(lightPosition[0][0], lightPosition[0][1]);//-2.0, 0.0); //light at the left
			break;
		case 2:
			model.setLightPosition(lightPosition[1][0], lightPosition[1][1]);//0.0, 2.0); //light at the top
			break;
		case 3:
			model.setLightPosition(lightPosition[2][0], lightPosition[2][1]);//0.0, -2.0); //light at the bottom
			break;
		case 4:
			model.setLightPosition(lightPosition[3][0], lightPosition[3][1]);//2.0, 0.0); //light at the right
			break;
		default:
			throw new RuntimeException("Error light number!");
		}
	}
	
	/**
	 * Solve the forward problem with given parameters
	 * @param lightNum
	 * @param _mu_a
	 * @param _mu_sp
	 * @param postfix
	 * @return
	 */
	public Vector solveForward(int lightNum, Function _mu_a, Function _mu_sp, String postfix) {
		setLightPosition(lightNum);
		model.setMu_a(_mu_a);
		model.setMu_sp(_mu_sp);
		
		//Solve Neumann problem on the big mesh
		Vector uBig = model.solveNeumann(meshBig);
		if(debug) Tools.plotVector(meshBig, outputFolder, "u_big__"+postfix+".dat", uBig);
		
		Vector uSmall = Tools.extractData(meshBig, meshSmall, uBig);
		if(debug) Tools.plotVector(meshSmall, outputFolder, "u_samll__"+postfix+".dat", uSmall);
		return uSmall;
	}
	
	/**
	 * Solve the problem without inclusion (Background with constant mu_a and mu_sp)
	 * @param lightNum
	 * @param postfix
	 * @return
	 */
	public Vector solveBk(int lightNum, String postfix) {
		return solveForward(lightNum, ModelParam.getBkMu_a(), C(init_mu_sp), "bk_"+postfix);
	}
	
	/**
	 * Solve the problem with inclusion specified by member variable mu_a and mu_sp
	 * @param lightNum
	 * @param postfix
	 * @return
	 */
	public Vector solveInclusion(int lightNum, String postfix) {
		if(debug) {
			Tools.plotFunction(meshSmall, outputFolder, "mu_a_real_"+postfix+".dat", mu_a);
			Tools.plotFunction(meshSmall, outputFolder, "mu_sp_real_"+postfix+".dat", mu_sp);
		}
		return solveForward(lightNum, this.mu_a, this.mu_sp, "inc_"+postfix);
	}
	
	public Vector getTailDiff(int lightNum, Vector uSmallBk, Vector uSmallInc, String postfix) {
		switch(lightNum) {
		case 1:
			return getTailDiffRight(uSmallBk,uSmallInc,postfix); //light at the left
		case 2:
			return getTailDiffBottom(uSmallBk,uSmallInc,postfix); //light at the top
		case 3:
			throw new RuntimeException("Error light number 3!");
		case 4:
			throw new RuntimeException("Error light number 4!");
		default:
			throw new RuntimeException("Error light number "+lightNum+"!");
		}
	}
	/**
	 * Extend the difference of light intensity from right boundary to all the domain
	 * 
	 * @param uSmallBk -- Only values on right are needed
	 * @param uSmallInc -- Only values on right are needed
	 * @param postfix
	 * @return
	 */
	public Vector getTailDiffRight(Vector uSmallBk, Vector uSmallInc, String postfix) {
		//diff = uSmallInc - uSmallBk (negative)
		Vector diff = uSmallBk.copy();
		diff.axpy(-1.0, uSmallInc); 
		
		Vector2Function fuDiff = new Vector2Function(diff, meshSmall, "x", "y");
		Vector tail = uSmallInc.copy();
		
		//Extend the difference of light intensity from right boundary to all the domain
		Variable v = new Variable();
		NodeList nodes = meshSmall.getNodeList();
		for(int i=1; i<=nodes.size(); ++i) {
			v.set("x", 1.0); //depends on the position of light source
			v.set("y", nodes.at(i).coord(2));
			tail.set(i, fuDiff.value(v));
		}
		
		if(debug) Tools.plotVector(meshSmall, outputFolder, "tail_diff_"+postfix+".dat", tail);
		return tail;
	}
	
	/**
	 * Extend the difference of light intensity from bottom boundary to all the domain
	 * 
	 * @param uSmallBk -- Only values on bottom are needed
	 * @param uSmallInc -- Only values on bottom are needed
	 * @param postfix
	 * @return
	 */
	public Vector getTailDiffBottom(Vector uSmallBk, Vector uSmallInc, String postfix) {
		//diff = uSmallInc - uSmallBk (negative)
		Vector diff = uSmallBk.copy();
		diff.axpy(-1.0, uSmallInc);
		
		Vector2Function fuDiff = new Vector2Function(diff, meshSmall, "x", "y");
		Vector tail = uSmallInc.copy();
		
		Variable v = new Variable();
		NodeList nodes = meshSmall.getNodeList();
		for(int i=1; i<=nodes.size(); ++i) {
			v.set("x", nodes.at(i).coord(1)); 
			v.set("y", -1.0);//depends on the position of light source
			tail.set(i, fuDiff.value(v));
		}
		if(debug) Tools.plotVector(meshSmall, outputFolder, "tail_diff_"+postfix+".dat", tail);
		return tail;
	}

	public void runAllTest(int testNum) {
		debug = false;
		//test1();
		//test2();
		//test3();
		
		debug = true;
//		for(int i=0; i<5; i++) {
//			final int a = i;
//			System.out.println(a);
//		}
	
		if(testNum == 1)
			test11();
		else if(testNum == 2)
			test21();
		else if(testNum == 3)
			test31();
		else if(testNum == 4)
			test41();
		
		
	}
	
	//mu_sp error: sin(x)*sin(y)
	public void test1() {
		for(int i=2; i<=4; ++i) {
			double max_real_mu_a = i/10.0;
			this.mu_a = ModelParam.getMu_a(0.0, 0.0, 0.5, //(x,y;r)
					max_real_mu_a, //maxMu_a
					1); //type
			double err_s = -0.15;
			double err_e = 0.15;
			double step = 0.01;
			int j=0;
			for(double err = err_s; err <=err_e; err+=step) {
				mu_sp_error = err;
				this.mu_sp = new AbstractFunction("x","y") {
					@Override
					public double value(Variable v) {
						double x = v.get("x");
						double y = v.get("y");
						return init_mu_sp+
								mu_sp_error*init_mu_sp*Math.sin(2*Math.PI*x)*
								Math.sin(2*Math.PI*y);
					}
				};
				Tools.plotFunction(meshSmall, outputFolder, 
						String.format("test1_mu_sp__%02d_%02d.dat",i,j), 
						this.mu_sp);
				j++;
				Vector reconstruted_mu_a = _run();
				
				System.out.println(String.format("Test1: Error=%1.2f, max(mu_a_real)=%1.3f,max(mu_a_recon)=%1.3f",
						err, max_real_mu_a, reconstruted_mu_a.normInf()));
			}
		}
	}
	
	//mu_sp error in the location of mu_a
	public void test2() {
		final double radius = 0.5;
		final double x0 = 0.0;
		final double y0 = 0.0;
		for(int i=2; i<=4; ++i) {
			double max_real_mu_a = i/10.0;
			this.mu_a = ModelParam.getMu_a(x0, y0, radius, //(x,y;r)
					max_real_mu_a, //maxMu_a
					1); //type
			double err_s = -0.15;
			double err_e = 0.15;
			double step = 0.01;
			int j=0;
			for(double err = err_s; err <=err_e; err+=step) {
				final double rand = Math.random()*2.0-1.0;
				mu_sp_error = err+step*0.7*rand;
				if(Math.abs(err) < 1e-8) mu_sp_error = 0.0;
				this.mu_sp = new AbstractFunction("x","y") {
					@Override
					public double value(Variable v) {
						double x = v.get("x");
						double y = v.get("y");
						double dx = x - x0;
						double dy = y - y0;
						if(Math.sqrt(dx*dx+dy*dy) < radius) 
							return init_mu_sp+mu_sp_error*init_mu_sp*Math.cos((Math.PI/2)*Math.sqrt(dx*dx+dy*dy)/radius); 
						else
							return init_mu_sp;
					}
				};
				Tools.plotFunction(meshSmall, outputFolder, 
						String.format("test2_mu_sp__%02d_%02d.dat",i,j), 
						this.mu_sp);
				j++;
				Vector reconstruted_mu_a = _run();
				
				System.out.println(String.format("Test2: Error=%1.3f, max(mu_a_real)=%1.3f,max(mu_a_recon)=%1.3f",
						mu_sp_error, max_real_mu_a, reconstruted_mu_a.normInf()));
			}
		}
	}	
	
	//mu_sp error is constant
	public void test3() {
		for(int i=2; i<=4; ++i) {
			double max_real_mu_a = i/10.0;
			this.mu_a = ModelParam.getMu_a(0.0, 0.0, 0.5, //(x,y;r)
					max_real_mu_a, //maxMu_a
					1); //type
			double err_s = -0.15;
			double err_e = 0.15;
			double step = 0.01;
			int j=0;
			for(double err = err_s; err <=err_e; err+=step) {
				mu_sp_error = err;
				this.mu_sp = new AbstractFunction("x","y") {
					@Override
					public double value(Variable v) {
						return init_mu_sp+mu_sp_error*init_mu_sp;
					}
				};
				Tools.plotFunction(meshSmall, outputFolder, 
						String.format("test3_mu_sp__%02d_%02d.dat",i,j), 
						this.mu_sp);
				j++;
				Vector reconstruted_mu_a = _run();
				
				System.out.println(String.format("Test3: Error=%1.2f, max(mu_a_real)=%1.3f,max(mu_a_recon)=%1.3f",
						err, max_real_mu_a, reconstruted_mu_a.normInf()));
			}
		}
	}
	
	/**
	 * Compute u0(x,y) with constant mu_a and mu_sp
	 * Compute u(x,y) with non-constant mu_sp
	 *  
	 * @return
	 */
	private Vector _run() {
		String postfix = "_right";
		int lightNum = 1;
		Vector uSmallBk = solveBk(lightNum, postfix);
		Vector uSmall = solveInclusion(lightNum, postfix);
		//k=1/(3*mu_sp)
		Vector sol_a = Tools.solveParamInverse(meshSmall, uSmall, C0, C(1.0/(3.0*init_mu_sp)), ModelParam.getBkMu_a());
		if(debug) Tools.plotVector(meshSmall, outputFolder, "mu_a_sol_real"+postfix+".dat", sol_a);
		Vector diffRight = getTailDiffRight(uSmallBk, uSmall, postfix);
		
		postfix = "_bottom";
		lightNum = 2;
		uSmallBk = solveBk(lightNum, postfix);
		uSmall = solveInclusion(lightNum, postfix);
		//k=1/(3*mu_sp)
		sol_a = Tools.solveParamInverse(meshSmall, uSmall, C0, C(1.0/(3.0*init_mu_sp)), ModelParam.getBkMu_a());
		if(debug) Tools.plotVector(meshSmall, outputFolder, "mu_a_sol_real_"+postfix+".dat", sol_a);
		Vector diffBottom = getTailDiffBottom(uSmallBk, uSmall, postfix);
		
		//tailDiff = diffRight + diffBottom
		Vector tailDiff = diffRight.copy();
		tailDiff.axpy(1.0, diffBottom);
		if(debug) Tools.plotVector(meshSmall, outputFolder, "tail_diff.dat", tailDiff);
		//Cut 92.7%
		double max = 0.927*FMath.abs(tailDiff).normInf();
		for(int i=1; i<=tailDiff.getDim(); ++i) {
			if(Math.abs(tailDiff.get(i)) < max)
				tailDiff.set(i,0.0);
			else {
				double t = tailDiff.get(i);
				tailDiff.set(i,t>0?t-max:t+max);
			}
		}
		if(debug) Tools.plotVector(meshSmall, outputFolder, "tail_diff_cut.dat", tailDiff);
		
		//tailDiff用在4个光源的重构中
		String[] postfixs = {"left","top","bottom","right"};
		int[] lightNums = {1, 2, 3, 4};
		Vector sol_a_tail_sum = tailDiff.copy();
		sol_a_tail_sum.setAll(0.0);
		for(int i=0; i<lightNums.length; ++i) {
			uSmallBk = solveBk(lightNums[i], postfixs[i]);
			//k=1/(3*mu_sp)
			Vector tail = tailDiff.copy().add(uSmallBk);
			Vector sol_a_tail = Tools.solveParamInverse(meshSmall, tail, C0, C(1.0/(3.0*init_mu_sp)), ModelParam.getBkMu_a());
			if(debug) Tools.plotVector(meshSmall, outputFolder, "mu_a_sol_tail"+postfixs[i]+".dat", sol_a_tail);
			sol_a_tail_sum.axpy(1.0, sol_a_tail);
		}
		sol_a_tail_sum.scale(1.0/lightNums.length);
		
		Vector sol_a_tail_smooth = Utils.gaussSmooth(meshSmall, sol_a_tail_sum, 1, 0.5);
		sol_a_tail_smooth = Utils.gaussSmooth(meshSmall, sol_a_tail_smooth, 1, 0.5);
		//cut negative
		for(int i=1; i<sol_a_tail_smooth.getDim(); ++i) {
			if(sol_a_tail_smooth.get(i) < ModelParam.getBkMu_a().value()) 
				sol_a_tail_smooth.set(i,ModelParam.getBkMu_a().value());
		}
		
		Vector uSmallMeasurement = solveInclusion(lightNums[0], postfixs[0]);
		ReconResult rlt = enhanceMu_a(sol_a_tail_smooth, uSmallMeasurement, this.init_mu_sp, lightNums[0]);
		if(debug) Tools.plotVector(meshSmall, outputFolder, "mu_a_sol_tail_enhanced.dat", rlt.reconMu_a);
		return rlt.reconMu_a;
	}
	
	/**
	 * Compute 2-norm on the border against the light source indcated by lightNum
	 * 
	 * @param u
	 * @param v
	 * @param lightNum
	 * @return
	 */
	public double norm2OnBorder(Vector u, Vector v, int lightNum) {
		NodeList nodes = meshSmall.getNodeList();
		double sum = 0.0;
		switch(lightNum) {
		case 1:
			for(int i=1; i<=nodes.size(); ++i) { //light at left
				Node n = nodes.at(i);
				if(Math.abs(n.coord(1)-1.0)<Constant.meshEps)
					sum += (u.get(i)-v.get(i))*(u.get(i)-v.get(i));
			}
			break;
		case 2:
			//break;
		case 3:
			//break;
		case 4:
			//break;
		default:
			throw new RuntimeException("Error light number!");
		}
		return Math.sqrt(sum);
	}
	public double computeNorm2OnBorder(Vector mu_a, double mult, Vector uSmallMeasurement, double mu_sp, 
			int lightNum, int cnt) {
		//Scale mu_a by mult
		Vector t = mu_a.copy();
		t.shift(-ModelParam.getBkMu_a().value());
		t.scale(mult);
		t.shift(ModelParam.getBkMu_a().value());
		
		Vector mu_a_big = Tools.extendData(meshSmall, meshBig, t, 
				ModelParam.getBkMu_a().value());
		model.setMu_a(new Vector2Function(mu_a_big, meshBig, "x","y"));
		//if(debug) Tools.plotFunction(meshBig, outputFolder, "mu_a_iter_"+cnt+".dat", model.getMu_a());
		model.setMu_sp(C(mu_sp));
		Vector uBig = model.solveNeumann(meshBig);
		Vector uSmall = Tools.extractData(meshBig, meshSmall, uBig);
		//if(debug) Tools.plotVector(meshSmall, outputFolder, "u_samll_iter_"+cnt+".dat", uSmall);
		return norm2OnBorder(uSmall, uSmallMeasurement, lightNum);
	}
	
	public static class ReconResult {
		public Vector reconMu_a;
		public double errorNormOnBoundary;
	}
	
	public ReconResult enhanceMu_a(Vector mu_a, Vector uSmallMeasurement, double mu_sp, int lightNum) {
		if(debug) System.out.print("enhanceMu_a ");
		setLightPosition(lightNum);

		//Find the increase point of the norm
		double norm = Double.MAX_VALUE;
		double mult = 1.0;
		int cnt = 0;
		//if(debug) System.out.println("Increase mult:");
		while(true) {
			//if(debug) System.out.println("mult="+mult);
			if(debug) System.out.print(".");
			double tmp = computeNorm2OnBorder(mu_a, mult, uSmallMeasurement, mu_sp, lightNum, cnt++);
			if(tmp < norm)
				norm = tmp;
			else 
				break;
			mult *= 2.0;
		}
		//Golden Section Search
		double a = 1.0, b = mult;
		double r = (Math.sqrt(5)-1.0)/2.0; //0.618...
		double eps = 0.1;
		double x = a + r*(b-a);
		double y = a + r*r*(b-a);
		//if(debug) System.out.println("Golden Section Search (enhanceMu_a): ");
		double u = computeNorm2OnBorder(mu_a, x, uSmallMeasurement, mu_sp, lightNum, cnt++);
		double v = computeNorm2OnBorder(mu_a, y, uSmallMeasurement, mu_sp, lightNum, cnt++);
		while(b-a > eps) {
			//if(debug) System.out.println(String.format("[%1.2f, %1.2f]: Norm=%1.5f", a,b,u));
			if(debug) System.out.print(".");

			if(u > v) {
				b = x;
				x = y;
				u = v;
				y = a + r*r*(b-a);
				v = computeNorm2OnBorder(mu_a, y, uSmallMeasurement, mu_sp, lightNum, cnt++);
			} else {
				a = y;
				y = x;
				v = u;
				x = a + r*(b-a);
				u = computeNorm2OnBorder(mu_a, x, uSmallMeasurement, mu_sp, lightNum, cnt++);
			}
		}
		Vector t = mu_a.copy();
		t.shift(-ModelParam.getBkMu_a().value());
		for(int i=1; i<t.getDim(); ++i)
			t.set(i, t.get(i)*a*(0.85+Math.random()/6.67));
		t.shift(ModelParam.getBkMu_a().value());
		ReconResult rlt = new ReconResult();
		rlt.reconMu_a = t;
		rlt.errorNormOnBoundary = u;
		if(debug) {
			System.out.println("max(recon_mu_a)="+rlt.reconMu_a.normInf()+
					" errorNormOnBoundary="+rlt.errorNormOnBoundary+
					" volumn(recon_mu_a)="+computeVolumn(meshSmall, rlt.reconMu_a));
			
		}
		
		return rlt;
	}
	
	/**
	 * Firstly, we simulate the real background u00(x,y) by solving the equation with non-constant mu_sp
	 * 
	 * Secondly, we use a model with constant mu_a and mu_sp to match u00(x,y) by calibrating the parameter 
	 * k^2=3*mu_a*mu_sp. The calibrated solution u0(x,y) will be obtained.
	 * 
	 * Thirdly, we compute tail by using the measured data on the boundary:
	 *  1) Compute the difference between the measured data and background u0(x,y)
	 *  2) Cut off unnecessary data
	 *  3) Obtain the tail by adding the difference to the background 
	 *  
	 * Finally, solve the inverse problem to obtain a(x,y)
	 * 
	 * 
	 * @param calibratedUBk -- Index is lightNum
	 * @param simUInc -- Index is lightNum
	 */
	private ReconResult _run2(Vector[] calibratedUBk, double mu_sp, Vector[] simUInc, double cutPercent) {
		/** Test code segment
		mu_sp_error = 0.15;
		Function mu_sp_bk = new AbstractFunction("x","y") {
			@Override
			public double value(Variable v) {
				double x = v.get("x");
				double y = v.get("y");
				return init_mu_sp+
						mu_sp_error*init_mu_sp*Math.sin(2*Math.PI*x)*
						Math.sin(2*Math.PI*y);
			}
		};
		int lightNum = 1; //left
		Vector simUbk = simulateBackground(lightNum,mu_sp_bk,"left");
		double targetRatio = findMaxMinRatio(lightNum, this.meshSmall, simUbk);
		System.out.println("targetRatio="+targetRatio);

		double findMu_sp = calibrateMu_sp(lightNum,targetRatio,"");
		System.out.println("Found mu_sp="+findMu_sp);
		*/
		
		/**************888 Test for log(I0/I)=k*d*mua*************
		//String[] lightPostfixs = {"left","top","bottom","right"};
		Vector quotient = FMath.log(simUInc[0].copy().axDivy(1.0, calibratedUBk[0]));
		Vector2Function fquotient = new Vector2Function(quotient, meshSmall, "x", "y");

		double lightX = -2.0;
		double lightY = 0.0;
		Vector mua = new SpaceVector(meshSmall.getNodeList().size());
		Vector mua_quasi = new SpaceVector(meshSmall.getNodeList().size());
		//Extend the difference of light intensity from right boundary to all the domain
		Variable v = new Variable();
		NodeList nodes = meshSmall.getNodeList();
		for(int i=1; i<=nodes.size(); ++i) {
			v.set("x", 1.0); //depends on the position of light source
			v.set("y", nodes.at(i).coord(2));
			double dx = lightX-v.get("x");
			double dy = lightY-v.get("y");
			double distance = Math.sqrt(dx*dx+dy*dy);
			double tmp = fquotient.value(v)/distance;
			mua.set(i, tmp);
		}
		if(debug) Tools.plotVector(meshSmall, outputFolder, "mua_light_left_.dat", mua);
		quotient = FMath.log(simUInc[1].copy().axDivy(1.0, calibratedUBk[1]));
		fquotient = new Vector2Function(quotient, meshSmall, "x", "y");
		lightX = 0.0;
		lightY = 2.0;
		for(int i=1; i<=nodes.size(); ++i) {
			v.set("x", nodes.at(i).coord(1)); //depends on the position of light source
			v.set("y", -0.5);
			double dx = lightX-v.get("x");
			double dy = lightY-v.get("y");
			double distance = Math.sqrt(dx*dx+dy*dy);
			double tmp = fquotient.value(v)/distance;
			mua.set(i, tmp);
			
			dx = lightX-nodes.at(i).coord(1);
			dy = lightY-nodes.at(i).coord(2);
			distance = Math.sqrt(dx*dx+dy*dy);
			if(distance > 0.3)
				mua_quasi.set(i, fquotient.value(v)/distance);
			else
				mua_quasi.set(i, 0.0);
		}
		if(debug) Tools.plotVector(meshSmall, outputFolder, "mua_light_top_.dat", mua);
		if(debug) Tools.plotVector(meshSmall, outputFolder, "mua_light_top_quasi_.dat", mua_quasi);
		*/
		
				
		
		Vector diffRight = getTailDiff(1, calibratedUBk[0], simUInc[0], "right");
		Vector diffBottom = getTailDiff(2, calibratedUBk[1], simUInc[1], "bottom");

		//tailDiff = diffRight + diffBottom
		Vector tailDiff = diffRight.copy();
		tailDiff.axpy(1.0, diffBottom);
		if(debug) Tools.plotVector(meshSmall, outputFolder, "tail_diff_avg.dat", tailDiff);
		//tailDiff should be negative
		tailDiff.shift(-tailDiff.normInf());
		if(debug) Tools.plotVector(meshSmall, outputFolder, "tail_diff_avg_shift.dat", tailDiff);
		//for(int i=1; i<=tailDiff.getDim(); ++i) {
		//	if(tailDiff.get(i) > 0.0) tailDiff.set(i, 0.0);
		//}
		//Cut 92.7%  0.938
		double max = cutPercent*FMath.abs(tailDiff).normInf();
		for(int i=1; i<=tailDiff.getDim(); ++i) {
			if(Math.abs(tailDiff.get(i)) < max)
				tailDiff.set(i,0.0);
			else {
				double t = tailDiff.get(i);
				tailDiff.set(i,t>0?t-max:t+max);
			}
		}
		if(debug) Tools.plotVector(meshSmall, outputFolder, "tail_diff_avg_cut.dat", tailDiff);

			
		//tailDiff用在4个光源的重构中，mu_a取4个的均值
		Vector sol_a_tail_sum = tailDiff.copy();
		sol_a_tail_sum.setAll(0.0);
		for(int i=0; i<lightNums.length; ++i) {
			Vector tail = tailDiff.copy().add(calibratedUBk[i]);
			//k=1/(3*mu_sp)
			Vector sol_a_tail = Tools.solveParamInverse(meshSmall, tail, C0, C(1.0/(3.0*mu_sp)), ModelParam.getBkMu_a());
			if(debug) Tools.plotVector(meshSmall, outputFolder, "mu_a_sol_tail_"+lightPostfixs[i]+".dat", sol_a_tail);
			sol_a_tail_sum.axpy(1.0, sol_a_tail);
		}
		sol_a_tail_sum.scale(1.0/lightNums.length);
		
		Vector sol_a_tail_smooth = Utils.gaussSmooth(meshSmall, sol_a_tail_sum, 1, 0.5);
		sol_a_tail_smooth = Utils.gaussSmooth(meshSmall, sol_a_tail_smooth, 1, 0.5);
		//cut negative
		for(int i=1; i<sol_a_tail_smooth.getDim(); ++i) {
			if(sol_a_tail_smooth.get(i) < ModelParam.getBkMu_a().value()) 
				sol_a_tail_smooth.set(i,ModelParam.getBkMu_a().value());
		}
		
		//Enhance mu_a （使用lightNum=1）
		ReconResult rlt = enhanceMu_a(sol_a_tail_smooth, simUInc[0], mu_sp, lightNums[0]);
		if(debug) Tools.plotVector(meshSmall, outputFolder, "mu_a_sol_tail_enhanced.dat", rlt.reconMu_a);
		return rlt;

	}
	
	//mu_s has a sin error on the whole domain
	private void test11() {
		System.out.println("Running test11()...");
		final double radius = 0.5;
		final double x0 = 0.0;
		final double y0 = 0.0;
		double err_s = -0.15;
		double err_e = 0.15;
		double step = 0.01;
		int j=0;
		for(double err = err_s; err <=err_e; err+=step) {
			mu_sp_error = err;
			this.mu_sp = new AbstractFunction("x","y") {
				@Override
				public double value(Variable v) {
					double x = v.get("x");
					double y = v.get("y");
					return init_mu_sp+
							mu_sp_error*init_mu_sp*Math.sin(2*Math.PI*x)*
							Math.sin(2*Math.PI*y);
				}
			};
			Tools.plotFunction(meshSmall, outputFolder, String.format("test11_mu_sp__%02d.dat",j++), this.mu_sp);
			System.out.println("Real extrema(mu_sp)="+(init_mu_sp+mu_sp_error*init_mu_sp));
			
			int lightNum = 1; //left
			String postfix = "left";
			Vector simUbk = simulateBackground(lightNum,mu_sp,postfix);
			double targetRatio = findMaxMinRatio(lightNum, this.meshSmall, simUbk);
			//System.out.println("Target Ratio="+targetRatio);

			final double calibratedMu_sp = calibrateMu_sp(lightNum,targetRatio,"");
			final Vector[] calibratedUBk = new Vector[4];
			for(int k=0; k<lightNums.length; ++k)
				calibratedUBk[k] = solveForward(lightNums[k], ModelParam.getBkMu_a(), C(calibratedMu_sp), lightPostfixs[k]);

			//mu_a = 0.2, 0.3, 0.4
			for(int i=2; i<=4; ++i) {
				double max_real_mu_a = i/10.0;
				this.mu_a = ModelParam.getMu_a(x0, y0, radius, //(x,y;r)
						max_real_mu_a, //maxMu_a
						1); //type
				System.out.println("volumn(real mu_a)="+computeVolumn(this.meshSmall, this.mu_a));
				
				final Vector[] simUInc = new Vector[4];
				for(int k=0; k<lightNums.length; ++k)
					simUInc[k] = solveInclusion(lightNums[k], lightPostfixs[k]);

				ReconResult rlt = _run2(calibratedUBk, calibratedMu_sp, simUInc, 0.90);//0.93
				System.out.println(String.format("Test11: Error=%1.3f, max(mu_a_real)=%1.3f,max(mu_a_recon)=%1.3f, volumn(mu_a)=%1.5f",
						mu_sp_error, 
						max_real_mu_a, 
						rlt.reconMu_a.normInf(),
						computeVolumn(this.meshSmall, rlt.reconMu_a))
						);
			}
		}		
	}
	
	//mu_s has an error at the location the same as the inclusion
	private void test21() {
		System.out.println("Running test21()...");
		final double radius = 0.5;
		final double x0 = 0.0;
		final double y0 = 0.0;
		double err_s = -0.15;
		double err_e = 0.15;
		double step = 0.01;
		int j=0;
		for(double err = err_s; err <=err_e; err+=step) {
			final double rand = Math.random()*2.0-1.0; //-1.0~1.0
			mu_sp_error = err+step*0.7*rand;
			if(Math.abs(err) < 1e-8) mu_sp_error = 0.0;
			this.mu_sp = new AbstractFunction("x","y") {
				@Override
				public double value(Variable v) {
					double x = v.get("x");
					double y = v.get("y");
					double dx = x - x0;
					double dy = y - y0;
					if(Math.sqrt(dx*dx+dy*dy) < radius) 
						return init_mu_sp+mu_sp_error*init_mu_sp*Math.cos((Math.PI/2)*Math.sqrt(dx*dx+dy*dy)/radius); 
					else
						return init_mu_sp;
				}
			};
			Tools.plotFunction(meshSmall, outputFolder, String.format("test21_mu_sp__%02d.dat",j++), this.mu_sp);
			System.out.println("Real extrema(mu_sp)="+(init_mu_sp+mu_sp_error*init_mu_sp));
			
			int lightNum = 1; //left
			String postfix = "left";
			Vector simUbk = simulateBackground(lightNum,mu_sp,postfix);
			double targetRatio = findMaxMinRatio(lightNum, this.meshSmall, simUbk);
			//System.out.println("Target Ratio="+targetRatio);

			final double calibratedMu_sp = calibrateMu_sp(lightNum,targetRatio,"");
			final Vector[] calibratedUBk = new Vector[4];
			for(int k=0; k<lightNums.length; ++k)
				calibratedUBk[k] = solveForward(lightNums[k], ModelParam.getBkMu_a(), C(calibratedMu_sp), lightPostfixs[k]);

			//mu_a = 0.2, 0.3, 0.4
			for(int i=2; i<=4; ++i) {
				double max_real_mu_a = i/10.0;
				this.mu_a = ModelParam.getMu_a(x0, y0, radius, //(x,y;r)
						max_real_mu_a, //maxMu_a
						1); //type
				System.out.println("volumn(real mu_a)="+computeVolumn(this.meshSmall, this.mu_a));
				
				final Vector[] simUInc = new Vector[4];
				for(int k=0; k<lightNums.length; ++k)
					simUInc[k] = solveInclusion(lightNums[k], lightPostfixs[k]);
				
//Golden Section Search for cut percentage
//				Function fun = new AbstractFunction("x") {
//					@Override
//					public double value(Variable v) {
//						double val = v.get();
//						ReconResult rlt = _run2(calibratedUBk, calibratedMu_sp, simUInc, val);
//						//return computeVolumn(meshSmall, rlt.reconMu_a);
//						return rlt.errorNormOnBoundary;
//					}
//				};
//				fun.setFName("_run2");
//				double cutPersent = Utils.GoldenSectionSearch(0.91, 0.98, 0.005, fun, true);
//				System.out.println("cut persent="+cutPersent);
//				ReconResult rlt = _run2(calibratedUBk, calibratedMu_sp, simUInc, cutPercent);
				
				ReconResult rlt = _run2(calibratedUBk, calibratedMu_sp, simUInc, 0.90);//0.93
				System.out.println(String.format("Test21: Error=%1.3f, max(mu_a_real)=%1.3f,max(mu_a_recon)=%1.3f, volumn(mu_a)=%1.5f",
						mu_sp_error, 
						max_real_mu_a, 
						rlt.reconMu_a.normInf(),
						computeVolumn(this.meshSmall, rlt.reconMu_a))
						);
			}
		}
	}
	
	//mu_s has a constant error in the whole domain
	//this will be discarded
	private void test31() {
		System.out.println("Running test31()...");
		final double radius = 0.5;
		final double x0 = 0.0;
		final double y0 = 0.0;
		double err_s = -0.05;
		double err_e = 0.05;
		double step = 0.01;
		int j=0;
		System.out.println("volumn(background mu_a)="+computeVolumn(this.meshSmall, ModelParam.getBkMu_a()));
		for(int i=2; i<=4; ++i) {
			double max_real_mu_a = i/10.0;
			Function mu_a  = ModelParam.getMu_a(x0, y0, radius, //(x,y;r)
					max_real_mu_a, //maxMu_a
					1); //type
			System.out.println("max(mu_a)="+max_real_mu_a+" volumn(real mu_a)="+computeVolumn(this.meshSmall, mu_a));
		}
		
		for(double err = err_s; err <=err_e; err+=step) {
			mu_sp_error = err;
			this.mu_sp = new AbstractFunction("x","y") {
				@Override
				public double value(Variable v) {
					return init_mu_sp+mu_sp_error*init_mu_sp;
				}
			};
			Tools.plotFunction(meshSmall, outputFolder, String.format("test31_mu_sp__%02d.dat",j++), this.mu_sp);
			System.out.println("Real extrema(mu_sp)="+(init_mu_sp+mu_sp_error*init_mu_sp));
			
			final Vector[] calibratedUBk = new Vector[4];
			for(int k=0; k<lightNums.length; ++k)
				calibratedUBk[k] = solveForward(lightNums[k], ModelParam.getBkMu_a(), C(init_mu_sp), lightPostfixs[k]);

			//mu_a = 0.2, 0.3, 0.4
			for(int i=2; i<=4; ++i) {
				double max_real_mu_a = i/10.0;
				this.mu_a = ModelParam.getMu_a(x0, y0, radius, //(x,y;r)
						max_real_mu_a, //maxMu_a
						1); //type
				System.out.println("volumn(real mu_a)="+computeVolumn(this.meshSmall, this.mu_a));
				
				final Vector[] simUInc = new Vector[4];
				for(int k=0; k<lightNums.length; ++k)
					simUInc[k] = solveInclusion(lightNums[k], lightPostfixs[k]);

				ReconResult rlt = _run2(calibratedUBk, init_mu_sp, simUInc, 0.90);//0.93
				System.out.println(String.format("Test31: Error=%1.3f, max(mu_a_real)=%1.3f,max(mu_a_recon)=%1.3f, volumn(mu_a)=%1.5f",
						mu_sp_error, 
						max_real_mu_a, 
						rlt.reconMu_a.normInf(),
						computeVolumn(this.meshSmall, rlt.reconMu_a))
						);
			}
		}		
	}
	
	
	//mu_s has a random error on the whole domain
	private void test41() {
		System.out.println("Running test41()...(random error on the whole domain)");
		final double radius = 0.5;
		final double x0 = 0.0;
		final double y0 = 0.0;
		double err_s = -0.15;
		double err_e = 0.15;
		double step = 0.01;
		int j=0;
		for(double err = err_s; err <=err_e; err+=step) {
			mu_sp_error = err;
			this.mu_sp = new AbstractFunction("x","y") {
				@Override
				public double value(Variable v) {
					return init_mu_sp+
							mu_sp_error*(init_mu_sp*1.2)*(Math.random()*2-1);
				}
			};
			Tools.plotFunction(meshSmall, outputFolder, String.format("test41_mu_sp__%02d.dat",j++), this.mu_sp);
			Vector mu_spv = Tools.function2vector(meshSmall, this.mu_sp);
			mu_spv = Utils.gaussSmooth(meshSmall, mu_spv, 1, 0.5);
			System.out.println("Real extrema(mu_sp)="+FMath.min(mu_spv));
			System.out.println("Real extrema(mu_sp)="+FMath.max(mu_spv));
			this.mu_sp = new Vector2Function(mu_spv, meshSmall, "x","y").setDefaultFunction(C(init_mu_sp));
			Tools.plotFunction(meshSmall, outputFolder, String.format("test41_mu_sp_smooth__%02d.dat",j++), this.mu_sp);
			
			int lightNum = 1; //left
			String postfix = "left";
			Vector simUbk = simulateBackground(lightNum,mu_sp,postfix);
			double targetRatio = findMaxMinRatio(lightNum, this.meshSmall, simUbk);
			//System.out.println("Target Ratio="+targetRatio);

//! test for log(I0/I)=k*d*mua
//!			final double calibratedMu_sp = this.init_mu_sp; 
			final double calibratedMu_sp = calibrateMu_sp(lightNum,targetRatio,"");
			final Vector[] calibratedUBk = new Vector[4];
			for(int k=0; k<lightNums.length; ++k)
				calibratedUBk[k] = solveForward(lightNums[k], ModelParam.getBkMu_a(), C(calibratedMu_sp), lightPostfixs[k]);

			//mu_a = 0.2, 0.3, 0.4
			for(int i=2; i<=4; ++i) {
//!			for(int i=3; i<=3; ++i) {
				double max_real_mu_a = i/10.0;
				this.mu_a = ModelParam.getMu_a(x0, y0, radius, //(x,y;r)
				max_real_mu_a, //maxMu_a
						1); //type
//!				this.mu_a = ModelParam.getMu_a(-0.25, 0.0, 0.3, //(x,y;r)
//!						max_real_mu_a, //maxMu_a
//!						2);//2, 0.5, false);
				System.out.println("volumn(real mu_a)="+computeVolumn(this.meshSmall, this.mu_a));
				
				final Vector[] simUInc = new Vector[4];
				for(int k=0; k<lightNums.length; ++k)
					simUInc[k] = solveInclusion(lightNums[k], lightPostfixs[k]);

				ReconResult rlt = _run2(calibratedUBk, calibratedMu_sp, simUInc, 0.90);//0.93
				System.out.println(String.format("Test41: Error=%1.3f, max(mu_a_real)=%1.3f,max(mu_a_recon)=%1.3f, volumn(mu_a)=%1.5f",
						mu_sp_error, 
						max_real_mu_a, 
						rlt.reconMu_a.normInf(),
						computeVolumn(this.meshSmall, rlt.reconMu_a))
						);
			}
		}
	}
	
	public double[] getCoord(int lightNum, Node node) {
//		String[] lightPostfixs = {"left","top","bottom","right"};
		if(lightNum==1) return new double[]{0.5, node.coord(2)};
		if(lightNum==2) return new double[]{node.coord(1), -0.3};
		if(lightNum==3) return new double[]{node.coord(1), 0.3};
		if(lightNum==4) return new double[]{-0.5, node.coord(2)};
		return null;
	}
	
	public void testBeerLambertLaw(int runCounter) {
		Vector[] mua_sol = new Vector[4];
		Vector uBigIncTop = null;
		for(int lightNum = 1; lightNum<=4; ++lightNum) {
			String postfix = this.lightPostfixs[lightNum-1];
			setLightPosition(lightNum);
			
			//Solve for background light intensity
			model.setMu_a(ModelParam.getBkMu_a());
			model.setMu_sp(C(this.init_mu_sp));
			//Solve Neumann problem on the big mesh
			Vector uBigBk = model.solveNeumann(meshBig);
			if(debug) Tools.plotVector(meshBig, outputFolder, "BeerLambertLaw_u_bigBk__"+postfix+".dat", uBigBk);
			
			//Solve for inclusion light intensity
			model.setMu_a(ModelParam.getMu_a(-0.25, 0.0, 0.3, //(x,y;r)
					0.3, //maxMu_a
					2));
			model.setMu_sp(C(this.init_mu_sp));
			mu_sp_error = 0.07; //err;
			model.setMu_sp(new AbstractFunction("x","y") {
				@Override
				public double value(Variable v) {
					return init_mu_sp+
							mu_sp_error*(init_mu_sp)*(Math.random()*2-1);
				}
			});
			//Solve Neumann problem on the big mesh
			Vector uBigInc = model.solveNeumann(meshBig);
			if(debug) Tools.plotVector(meshBig, outputFolder, "BeerLambertLaw_u_bigInc__"+postfix+".dat", uBigInc);
			if(lightNum == 1) uBigIncTop = uBigInc;
			
			Vector quotient = FMath.log10(uBigBk.copy().axDivy(1.0, uBigInc));
	
			double lightX = this.lightPosition[lightNum-1][0];
			double lightY = this.lightPosition[lightNum-1][1];
			Vector mua_quasi = new SpaceVector(meshBig.getNodeList().size());
			NodeList nodes = meshBig.getNodeList();
			for(int i=1; i<=nodes.size(); ++i) {
				double dx = lightX-nodes.at(i).coord(1);
				double dy = lightY-nodes.at(i).coord(2);
				double distance = Math.sqrt(dx*dx+dy*dy);
				if(distance > 0.3)
					mua_quasi.set(i, quotient.get(i)/distance);
				else
					mua_quasi.set(i, 0.0);
			}
			if(debug) Tools.plotVector(meshBig, outputFolder, 
					"BeerLambertLaw_mua_quasi_light_"+postfix+"_.dat", mua_quasi);

			//
			Vector mua_quasi_small = Tools.extractData(meshBig, meshSmall, mua_quasi);
			if(debug) Tools.plotVector(meshSmall, outputFolder, 
					"BeerLambertLaw_mua_quasi_small_light_"+postfix+".dat", mua_quasi_small);
			Vector mua = mua_quasi_small.copy();
			Vector mua_times_distance = mua_quasi_small.copy();
			mua.setAll(0.0);
			mua_times_distance.setAll(0.0);
			Vector2Function fmua_quasi_small = new Vector2Function(mua_quasi_small, meshSmall, "x", "y");
			NodeList nodesSmall = meshSmall.getNodeList();
			Variable v = new Variable();
			for(int i=1; i<=nodesSmall.size(); ++i) {
				//depends on the position of light source
				double[] coord = getCoord(lightNum, nodesSmall.at(i));
				v.set("x", coord[0]); 
				v.set("y", coord[1]);
				double dx = lightX-nodesSmall.at(i).coord(1);
				double dy = lightY-nodesSmall.at(i).coord(2);
				double distance = Math.sqrt(dx*dx+dy*dy);
				
				double tmp = fmua_quasi_small.value(v);
				mua.set(i, tmp);
				mua_times_distance.set(i, tmp*distance);
			}
			if(debug) Tools.plotVector(meshSmall, outputFolder, 
					"BeerLambertLaw_mua_recon_light_"+postfix+".dat", mua);
			
			//log(I0/I)=k*d*mua
			Vector uSmallBk = Tools.extractData(meshBig, meshSmall, uBigBk);
			
			Vector uGauss = uSmallBk.copy().axDivy(1.0, FMath.pow(10.0, mua_times_distance));
			mua_sol[lightNum-1] = Tools.solveParamInverse(meshSmall, uGauss, C0, C(1.0/(3.0*this.init_mu_sp)), ModelParam.getBkMu_a());
			if(debug) Tools.plotVector(meshSmall, outputFolder, 
					"BeerLambertLaw_mua_sol_light_"+postfix+".dat", mua_sol[lightNum-1]);
		}
		
		
		//tailDiff用在4个光源的重构中，mu_a取4个的均值
		Vector sol_a_tail_sum = mua_sol[0].copy();
		sol_a_tail_sum.setAll(0.0);
		Tools.scaleVectorsByFirst(mua_sol[0], mua_sol);
		for(int i=0; i<lightNums.length; ++i) {
			sol_a_tail_sum.axpy(1.0, mua_sol[i]);
		}
		sol_a_tail_sum.scale(1.0/lightNums.length);
		
		Vector sol_a_tail_smooth = Utils.gaussSmooth(meshSmall, sol_a_tail_sum, 2, 0.5);
		sol_a_tail_smooth = Utils.gaussSmooth(meshSmall, sol_a_tail_smooth, 2, 0.5);
		//cut negative
		for(int i=1; i<sol_a_tail_smooth.getDim(); ++i) {
			if(sol_a_tail_smooth.get(i) < ModelParam.getBkMu_a().value()) 
				sol_a_tail_smooth.set(i,ModelParam.getBkMu_a().value());
		}
		if(debug) Tools.plotVector(meshSmall, outputFolder, "BeerLambertLaw_mu_a_sol_all.dat", 
				sol_a_tail_smooth);
		//cut
		double max = 0.827*FMath.max(sol_a_tail_smooth);
		for(int i=1; i<=sol_a_tail_smooth.getDim(); ++i) {
			if(sol_a_tail_smooth.get(i) < max)
				sol_a_tail_smooth.set(i, 0.1);//background mua
		}
		if(debug) Tools.plotVector(meshSmall, outputFolder, "BeerLambertLaw_mu_a_sol_all_cut.dat", 
				sol_a_tail_smooth);
		
		//Enhance mu_a （使用lightNum=1）
		Vector uSmallInc = Tools.extractData(meshBig, meshSmall, uBigIncTop);
		ReconResult rlt = enhanceMu_a(sol_a_tail_smooth,uSmallInc, init_mu_sp, lightNums[0]);
		Tools.plotVector(meshSmall, outputFolder, "BeerLambertLaw_mu_a_sol_all_cut_enhanced"+runCounter+".dat", rlt.reconMu_a);

	}
	
	public double computeVolumn(Mesh mesh, Function f) {
		double sum = 0.0;
		NodeList nodes = mesh.getNodeList();
		Variable v = new Variable();
		for(int i=1; i<nodes.size(); ++i) {
			Node n = nodes.at(i);
			v.set("x", n.coord(1));
			v.set("y", n.coord(2));
			sum += f.value(v);
		}
		double elementArea = 0.1*0.1;
		return sum*elementArea;
	}
	
	public double computeVolumn(Mesh mesh, Vector v) {
		double sum = 0.0;
		for(int i=1; i<v.getDim(); ++i) {
			sum += v.get(i);
		}
		double elementArea = 0.1*0.1;
		return sum*elementArea;
	}
	
	/**
	 * 
	 * @param lightNum
	 * @param postfix
	 * @return
	 */
	public Vector simulateBackground(int lightNum, Function mu_sp_bk, String postfix) {
		setLightPosition(lightNum);
		//--------------Problem without inclusion------------------------------------
		model.setMu_a(ModelParam.getBkMu_a());
		if(debug) Tools.plotFunction(meshSmall, outputFolder, "mu_a_bk_"+postfix+".dat", model.getMu_a());
		model.setMu_sp(mu_sp_bk);
		if(debug) Tools.plotFunction(meshSmall, outputFolder, "mu_sp_bk_"+postfix+".dat", model.getMu_sp());
		
		//Solve neumann problem on big mesh
		Vector uBig = model.solveNeumann(meshBig);
		if(debug) Tools.plotVector(meshBig, outputFolder, "u_big_bk_"+postfix+".dat", uBig);
		
		Vector uSmall = Tools.extractData(meshBig, meshSmall, uBig);
		if(debug) Tools.plotVector(meshSmall, outputFolder, "u_samll_bk_"+postfix+".dat", uSmall);
		return uSmall;
	}
	
	
	/**
	 * Return calibrated mu_sp
	 * @param lightNum
	 * @param maxMinRatio
	 * @return
	 */
	public double calibrateMu_sp(int lightNum, double maxMinRatio, String postfix) {
		if(debug) System.out.print("calibrateMu_sp ");
		//Golden Section Search
		int cnt = 0;
		double a = init_mu_sp/2, b = init_mu_sp*2;
		double r = (Math.sqrt(5)-1.0)/2.0; //0.618...
		double eps = 0.01;
		double x = a + r*(b-a);
		double y = a + r*r*(b-a);
		//if(debug) System.out.println("Golden Section Search (calibrateMu_sp):");
		double u = Math.abs(computeBackgroundMaxMinRatio(lightNum, x, "iter"+(cnt++)+"")-maxMinRatio);
		double v = Math.abs(computeBackgroundMaxMinRatio(lightNum, y, "iter"+(cnt++)+"")-maxMinRatio);
		while(b-a > eps) {
			//if(debug) System.out.println(String.format("[%1.2f, %1.2f]: Abs(targetRatio-currentRatio)=%1.5f", a,b,u));
			if(debug) System.out.print(".");
			if(u > v) {
				b = x;
				x = y;
				u = v;
				y = a + r*r*(b-a);
				v = Math.abs(computeBackgroundMaxMinRatio(lightNum, y, "iter"+(cnt++)+"")-maxMinRatio);
			} else {
				a = y;
				y = x;
				v = u;
				x = a + r*(b-a);
				u = Math.abs(computeBackgroundMaxMinRatio(lightNum, x, "iter"+(cnt++)+"")-maxMinRatio);
			}
		}
		System.out.println("Calibrated mu_sp="+x);
		return x;
	}
	
	public double computeBackgroundMaxMinRatio(int lightNum, double mu_sp, String postfix) {
		setLightPosition(lightNum);
		model.setMu_a(ModelParam.getBkMu_a());
		model.setMu_sp(C(mu_sp));
		
		//Solve neumann problem on big mesh
		Vector uBigBk = model.solveNeumann(meshBig);
		if(debug) Tools.plotVector(meshBig, outputFolder, "u_big_bk_"+postfix+".dat", uBigBk);
		
		Vector uSmallBk = Tools.extractData(meshBig, meshSmall, uBigBk);
		if(debug) Tools.plotVector(meshSmall, outputFolder, "u_samll_bk_"+postfix+".dat", uSmallBk);
		
		return findMaxMinRatio(lightNum, this.meshSmall, uSmallBk);
	}	
	
	
	/**
	 * Find out the ratio between the maximum and the minimum 
	 * value of the light intensity u on the given mesh.
	 */
	public double findMaxMinRatio(int lightNum, Mesh mesh, Vector u) {
		double[] coordMax = {0,0}, coordMin={0,0};
		if(lightNum == 1) { //left
			coordMax[0] = -1; coordMax[1] = 0;
			coordMin[0] =  1; coordMin[1] = 0;
		} else if(lightNum == 2) { //top
			coordMax[0] = 0; coordMax[1] =  1;
			coordMin[0] = 0; coordMin[1] = -1;
		} else if(lightNum == 3) { //bottom
			coordMax[0] = 0; coordMax[1] = -1;
			coordMin[0] = 0; coordMin[1] =  1;
		} else if(lightNum == 4) { //right
			coordMax[0] =  1; coordMax[1] = 0;
			coordMin[0] = -1; coordMin[1] = 0;
		}
		Node nMax = mesh.findNodeByCoord(coordMax, 0.05);
		Node nMin = mesh.findNodeByCoord(coordMin, 0.05);
		double ratio = 
				computeMeanValueNearANode(nMax,u)/computeMeanValueNearANode(nMin,u);
		return ratio;
	}
	
	public double computeMeanValueNearANode(Node node, Vector u) {
		NodeList neibors = node.neighbors;
		double sum = u.get(node.getIndex());
		for(Node n : neibors) {
			sum += u.get(n.getIndex());
		}
		return sum/(neibors.size()+1);
	}
	
	/**
	 * Get forward matrix A
	 * Let the absorption coefficient has a small change at a node each time,
	 * then the following values on the boundary:
	 *     log(light intensity/background light intensity)
	 * form a column of A.
	 * 
	 * @return
	 */
	public Matrix getForwardMatrixA() {
		//Number of measurement for one light source
		int nMea = getNumberOfMeasurement();
		//Dimension of the forward matrix
		int nRows = nMea*lightNums.length;
		int nCols = meshSmall.getNodeList().size();
		
		//TODO Forward matrix is a full matrix, but class FullMatrix is not suitable
		//for the purpose (是否考虑增加类似matlab中的full matrix Java类？)
		Matrix A = new SparseMatrixRowMajor(nRows, nCols);
		for(int i=0; i<lightNums.length; ++i) {
			System.out.println("Processing light number: "+lightNums[i]);
			System.out.println(nCols+" columns:");
			for(int j=1; j<=nCols; ++j) {
			//for(int j=1; j<=5; ++j) {
				Vector v = getColumnOfA(lightNums[i], j);
				for(int k=1; k<=v.getDim(); ++k) {
					A.set(i*nMea+k, j, v.get(k));
				}
				System.out.print(j+" ");
			}
			System.out.println();
		}
		MatlabMatFileWriter w = new MatlabMatFileWriter();
		w.addMatrix(A);
		w.writeFile(this.outputFolder+"/forwardA2.mat");
		System.out.println("Write forward matrix A done!");
		return A;
	}
	
	public Matrix getForwardMatrixAFast() {
		int nNodes = meshSmall.getNodeList().size();
		int nMea = getNumberOfMeasurement();
		//Dimension of the forward matrix
		int nRows = nMea*lightNums.length;
		int nCols = nNodes;
		
		Matrix A = new SparseMatrixRowMajor(nRows, nCols);
		for(int lightNum=0; lightNum<lightNums.length; ++lightNum) {
			System.out.println("Processing light number: "+lightNums[lightNum]);
			System.out.println(nCols+" columns:");
			double sz = 0.26;
			for(double y=1.0; y>-1.0; y-=sz) {
				for(double x=-1.0; x<1.0; x+=sz) {
					NodeList blockNodes = findNodes(meshBig, x, y, x+sz, y-sz);
					System.out.println("("+x+","+y+")-("+(x+sz)+","+(y-sz)+"): "+blockNodes);
					Vector v = getColumnOfAByBlock(lightNums[lightNum], blockNodes, "");
					//Set columns of A
					for(int j=1; j<blockNodes.size(); ++j) {
						int colIndex = meshSmall.findNode(blockNodes.at(j)).getIndex();
						for(int i=1; i<=v.getDim(); ++i) {
							int rowIndex = lightNum*nMea+i;
							A.set(rowIndex, colIndex, v.get(i));
						}
					}
				}
			}
			System.out.println();
		}
		MatlabMatFileWriter w = new MatlabMatFileWriter();
		w.addMatrix(A);
		w.writeFile(this.outputFolder+"/forwardA0.26.mat");
		System.out.println("Write forward matrix A done!");
		return A;
	}
	/**
	 * Return a list of nodes of the mesh which are located in the rectangle
	 * (x1,y1)-(x2,y2)
	 * @param mesh
	 * @param x1
	 * @param y1
	 * @param x2
	 * @param y2
	 * @return
	 */
	public NodeList findNodes(Mesh mesh, double x1, double y1, double x2, double y2) {
		NodeList nodes = mesh.getNodeList();
		NodeList rlt = new NodeList();
		for(int i=1; i<nodes.size(); ++i) {
			Node n = nodes.at(i);
			double x = n.coord(1);
			double y = n.coord(2);
			if(x1 <= x && x <= x2 && y2 <= y && y <= y1) {
				rlt.add(n);
			}
		}
		return rlt;
	}
	
	private HashMap<Integer, Vector> uBoundaryBk = new HashMap<Integer, Vector>();
	public Vector getColumnOfA(int lightNum, int column) {
		setLightPosition(lightNum);
		//Compute background if necessary
		if(uBoundaryBk.get(lightNum) == null) {
			model.setMu_a(ModelParam.getBkMu_a());
			model.setMu_sp(C(this.init_mu_sp));
			//Solve Neumann problem on the big mesh
			Vector uBigBk = model.solveNeumann(meshBig);
			uBoundaryBk.put(lightNum, extractBorderValues(meshBig, uBigBk));
		}
		
		//Set a small change for the absorption coefficient
		Vector vMu_a = new SparseVectorHashMap(meshBig.getNodeList().size());
		Node n = meshSmall.getNodeList().at(column);
		vMu_a.set(meshBig.findNode(n).getIndex(), 1);
		
		model.setMu_a(new Vector2Function(vMu_a, meshBig, "x","y"));
		if(debug) Tools.plotVector(meshBig, outputFolder, "Mua_light_"+lightNum+"_col"+column+".dat", vMu_a);
		model.setMu_sp(C(this.init_mu_sp));
		//Solve Neumann problem on the big mesh
		Vector uBig = model.solveNeumann(meshBig);
		if(debug) Tools.plotVector(meshBig, outputFolder, "uBig_light_"+lightNum+"_col"+column+".dat", uBig);
		
		Vector uBoundary = extractBorderValues(meshBig, uBig);
		return FMath.log(uBoundary.axDivy(1.0, uBoundaryBk.get(lightNum)));
	}
	
	/**
	 * 
	 * @param lightNum
	 * @param blockNodes
	 * @return
	 */
	public Vector getColumnOfAByBlock(int lightNum, NodeList blockNodes, String postfix) {
		setLightPosition(lightNum);
		if(uBoundaryBk.get(lightNum) == null) {
			model.setMu_a(ModelParam.getBkMu_a());
			model.setMu_sp(C(this.init_mu_sp));
			//Solve Neumann problem on the big mesh
			Vector uBigBk = model.solveNeumann(meshBig);
			uBoundaryBk.put(lightNum, extractBorderValues(meshBig, uBigBk));
		}
		
		//Set a small change for the absorption coefficient
		Vector vMu_a = new SparseVectorHashMap(meshBig.getNodeList().size());
		for(int i=1; i<=blockNodes.size(); ++i)
			vMu_a.set(blockNodes.at(i).getIndex(), 1);
		
		model.setMu_a(new Vector2Function(vMu_a, meshBig, "x","y"));
		if(debug) Tools.plotVector(meshBig, outputFolder, "Mua_light_"+lightNum+"_block"+postfix+".dat", vMu_a);
		model.setMu_sp(C(this.init_mu_sp));
		//Solve Neumann problem on the big mesh
		Vector uBig = model.solveNeumann(meshBig);
		if(debug) Tools.plotVector(meshBig, outputFolder, "uBig_light_"+lightNum+"_block"+postfix+".dat", uBig);
		
		Vector uBoundary = extractBorderValues(meshBig, uBig);
		return FMath.log(uBoundary.axDivy(1.0, uBoundaryBk.get(lightNum))).ax(1.0/blockNodes.size());
	}
	
	public Vector getMeasurementVector() {
		int nMea = getNumberOfMeasurement();
		double err_s = -0.7;
		double err_e = 0.7;
		double step = 0.1;
		int counter = 0;
		Vector vRlt = new SpaceVector(nMea*lightNums.length);
		for(double err = err_s; err <=err_e; err+=step) {
			mu_sp_error = 0.0; //0.3; //err;
			this.mu_sp = new AbstractFunction("x","y") {
				@Override
				public double value(Variable v) {
					return init_mu_sp+
							mu_sp_error*(init_mu_sp*1.2)*(Math.random()*2-1);
				}
			};
			Tools.plotFunction(meshSmall, outputFolder, String.format("Ticknov_mu_sp__%02d.dat",counter), this.mu_sp);
			Vector mu_spv = Tools.function2vector(meshSmall, this.mu_sp);
			mu_spv = Utils.gaussSmooth(meshSmall, mu_spv, 1, 0.5);
			System.out.println("Real extrema(mu_sp)="+FMath.min(mu_spv));
			System.out.println("Real extrema(mu_sp)="+FMath.max(mu_spv));
			this.mu_sp = new Vector2Function(mu_spv, meshSmall, "x","y").setDefaultFunction(C(init_mu_sp));
			Tools.plotFunction(meshSmall, outputFolder, String.format("Ticknov_mu_sp_smooth__%02d.dat",counter), this.mu_sp);

			for(int i=0; i<lightNums.length; ++i) {
				//light position
				int lightNum = lightNums[i];
				this.setLightPosition(lightNum);
				
				//Background
				model.setMu_a(ModelParam.getBkMu_a());
				model.setMu_sp(C(this.init_mu_sp));
				//Solve Neumann problem on the big mesh
				Vector uBigBk = model.solveNeumann(meshBig);
				Vector uBigBoundaryBk = extractBorderValues(meshBig, uBigBk);
				
				//Inclusion
				this.mu_a = ModelParam.getMu_a(-0.25, 0, 0.3, //(x,y;r)
						0.3, //maxMu_a
						2); //type
				// this.mu_sp = C(this.init_mu_sp);
				Vector uSmall = solveInclusion(lightNum, this.lightPostfixs[i]);
				Vector uBigBoundaryMea = extractBorderValues(meshSmall, uSmall);
				
				Vector dis = getBorderDistanceVector(meshSmall, lightNum);
				Vector vDelta = FMath.log(uBigBoundaryMea.axDivy(1.0, uBigBoundaryBk));
				vDelta.axDivy(1.0, dis);
				
				for(int k=1; k<=nMea; ++k) {
					vRlt.set(i*nMea+k, vDelta.get(k));
				}
			}
			MatlabMatFileWriter w = new MatlabMatFileWriter();
			w.addVector(vRlt);
			w.writeFile(this.outputFolder+"/y"+counter+".mat");

			counter++;
		}
		System.out.println("Write measurement vector done!");
		return vRlt;
	}
	
	public void readReconstruction(String matFile, String varName) {
		MatlabMatFileReader r = new MatlabMatFileReader(this.outputFolder+"/"+matFile);
		Matrix m = r.getMatrix(varName);
		//See CVX.pdf
		//x_qp = quadprog( 2*A’*A, -2*A’*b, [], [], [], [], l, u );
		//Matrix m = r.getMatrix("x_qp");
		//Matrix m = r.getMatrix("x_lsq");
		Vector v = new SpaceVector(m.getRowDim());
		for(int i=1; i<=m.getColDim(); ++i) {
			for(int j=1; j<=m.getRowDim(); ++j)
				v.set(j,m.get(j, i));
			Tools.plotVector(meshSmall, outputFolder, "recon_mua_"+matFile+"_lambda=1e-"+i+".dat", v);
		}
	}
	
	
	/**
	 * Return the number of nodes on the boundary of meshSmall
	 * 
	 * @return
	 */
	public int getNumberOfMeasurement() {
		//4*21-4 = 80
		return meshSmall.getBoundaryNodeList().size();
	}
	
	public ArrayList<Node> getBorderNodes(Mesh mesh) {
		ArrayList<Node> list = new ArrayList<Node>();
		//top -> right -> bottom -> left
		for(int i=-10; i<=10; ++i) {
			Node n = mesh.findNode(new Node(0, i/10, 1.0));
			list.add(n);
		}
		for(int i=9; i>=-10; --i) {
			Node n = mesh.findNode(new Node(0, 1.0, i/10));
			list.add(n);
		}
		for(int i=9; i>=-10; --i) {
			Node n = mesh.findNode(new Node(0, i/10, -1.0));
			list.add(n);
		}
		for(int i=-9; i<=9; ++i) {
			Node n = mesh.findNode(new Node(0, -1.0, i/10));
			list.add(n);
		}
		return list;
	}
	
	public Vector getBorderDistanceVector(Mesh mesh, int lightNum) {
		ArrayList<Node> list = getBorderNodes(mesh);
		Vector v = new SpaceVector(list.size());
		for(int i=0; i<list.size(); ++i) {
			double lightX = this.lightPosition[lightNum-1][0];
			double lightY = this.lightPosition[lightNum-1][1];
			double dx = lightX - list.get(i).coord(1);
			double dy = lightY - list.get(i).coord(2);
			double distance = Math.sqrt(dx*dx+dy*dy);
			v.set(i+1, distance);
		}
		return v;		
	}
	
	/**
	 * Extract the data on the boundary of the domain [-1,1]*[-1,1]
	 * 
	 * @param mesh Can be this.meshBig or this.meshSmall
	 * @param u
	 * @return
	 */
	private Vector extractBorderValues(Mesh mesh, Vector u) {
		ArrayList<Node> list = getBorderNodes(mesh);
		Vector v = new SpaceVector(list.size());
		for(int i=0; i<list.size(); ++i) {
			v.set(i+1, u.get(list.get(i).getIndex()));
		}
		return v;
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		System.out.println("Running...");
		
		BOEPaper ins = new BOEPaper();
		//ins.getForwardMatrixA();
		//ins.getForwardMatrixAFast();
		//ins.getMeasurementVector();
		//for(int i=0; i<=0; ++i) ins.readReconstruction("X"+i+".mat", "X");
		
		//ins.runAllTest(4);
		
		ins.debug = false;
		for(int i=20; i<=40; ++i)
			ins.testBeerLambertLaw(i);
		
		/**
		if(args.length==0)
			System.out.println("No test number");
		else {
			if(args[0].equals("11"))
				ins.runAllTest(1);
			else if(args[0].equals("21"))
				ins.runAllTest(2);
			else if(args[0].equals("31"))
				ins.runAllTest(3);
			else if(args[0].equals("41"))
				ins.runAllTest(4);
			else
				System.out.println("Unknown param: "+args[0]);
		}
		*/
		System.out.println("Done!");
		StringBuilder sb;
		sb.de
		
	}

}
