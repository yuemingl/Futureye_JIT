package edu.uta.futureye.application;

import java.util.ArrayList;

import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.function.AbstractMathFunc;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.Vector2Function;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.io.MeshReader;
import static edu.uta.futureye.function.operator.FMath.*;

public class Curvature {
	String outputFolder = "Curvature";
	String mesh_file = "prostate_test9_ex.grd";
	double x0 = 2.0;
	double y0 = 2.2;
	double radius = 0.5;
	double mua = 0.2;
	int num_inclusions = 1;
	double distance = 0.0;
	
	Mesh meshBig = null;

	public void loadMesh() {
		//prostate_test9_ex.grd比prostate_test8_ex.grd更密一些
//		String gridFileBig = "prostate_test8_ex.grd";
		String gridFileBig = mesh_file;
		
		MeshReader readerForward = new MeshReader(gridFileBig);
		meshBig = readerForward.read2DMesh();
		
		//Use element library to assign degree of freedom (DOF) to element
		Tools.assignLinearShapFunction(meshBig);
		meshBig.computeNodeBelongsToElements();
		meshBig.computeNeighborNodes();
	}
	
	public static double[] range(double s, double e, int n) {
		double[] ret = new double[n];
		assert(n > 0);
		double delta = (e-s)/(n-1);
		ret[0] = s;
		for(int i=1; i<n; ++i) {
			ret[i] = ret[i-1] + delta;
		}
		return ret;
	}
	public void reverse(double[] ary) {
		int i=0, j=ary.length-1;
		while(i < j) {
			double t = ary[i];
			ary[i] = ary[j];
			ary[j] = t;
			i++;
			j--;
		}
	}
	
	public void run() {
		loadMesh();
	
		ModelDOT model = new ModelDOT();
		
		//double[] lx = {1,2,3,4,          4.5,4.5,4.5,      4,3,2,1,             0.5,0.5,0.5};
		//double[] ly = {3.5,3.5,3.5,3.5,  3,2,1,     0.5,0.5,0.5,0.5,0.5,  1,2,3};
		//Light positions form a rectangle
		ArrayList<Variable> list = new ArrayList<Variable>();
		double[] xr = range(0.5, 4.5, 17);
		double[] yr = range(0.5, 3.5, 13);
		for(double x : xr) list.add(new Variable("x",  x).set("y",3.5));//top
		reverse(yr);
		for(double y : yr) list.add(new Variable("x",4.5).set("y",  y));//right
		reverse(xr);
		for(double x : xr) list.add(new Variable("x",  x).set("y",0.5));//bottom
		reverse(yr);
		for(double y : yr) list.add(new Variable("x",0.5).set("y",  y));//left
		
		//for(int i=0;i<7;++i) list.remove(0);
		
		System.out.println("xy = [");
		for(int i=0; i<list.size(); ++i) {
			Variable xy = list.get(i);
			model.setLightPosition(xy.get("x"), xy.get("y"));

			//model.k = C1;
			model.setMu_a(ModelParam.getMu_a(x0, y0, radius, //(x,y;r)
					mua, //maxMu_a
					num_inclusions, distance, true)); //distance between two inclusions
			if(i==0)
				Tools.plotFunction(meshBig, outputFolder, "a_real.dat", model.getMu_a());
			Vector uBig = model.solveNeumann(meshBig);
			Tools.plotVector(meshBig, outputFolder, String.format("u_big_%02d.dat", i), uBig);
			
			model.setMu_a(C(0.1));
			Vector uBig0 = model.solveNeumann(meshBig);
			Tools.plotVector(meshBig, outputFolder, String.format("u_big0_%02d.dat", i), uBig0);
			
			Vector uDiff = uBig0.axpy(-1.0, uBig);
			Tools.plotVector(meshBig, outputFolder, String.format("u_diff_%02d.dat", i), uDiff);
			
			BoundaryValues fb = new BoundaryValues(uDiff);
			Tools.plotFunction(meshBig, outputFolder, String.format("u_diff_brdy_%02d.dat", i), fb);
			
			Variable pos = findPosition(uDiff);
			System.out.println(pos.get("x")+","+pos.get("y"));
		}
		System.out.println("];");
		System.out.println("scatter(xy(:,1),xy(:,2))");
	}
	
	//extract bounday values [1.0, 4.0]*[1.0, 3.0]
	class BoundaryValues extends AbstractMathFunc {
		MathFunc f = null;
		public BoundaryValues(Vector uDiff) {
			super("x","y");
			f = new Vector2Function(uDiff, meshBig, "x", "y");
		}
		@Override
		public double apply(Variable v) {
			double x = v.get("x");
			double y = v.get("y");
			double eps = 0.05;
			if(Math.abs(x-1.0)<eps || 
					   Math.abs(x-4.0)<eps || 
					   Math.abs(y-1.0)<eps || 
					   Math.abs(y-3.0)<eps)
				return f.apply(v);
			return 0;
		}
	}
	
	public Variable findPosition(Vector  uDiff) {
		MathFunc f = new Vector2Function(uDiff, meshBig, "x", "y");
		Variable ret = new Variable();
		double delta = 0.02;
		Variable v = new Variable();

		//Find out the min value on top boundary
		double min = Double.MAX_VALUE;
		double x = 0.0;
		v.set("y",3.0);
		while(x < 5.0) {
			v.set("x", x);
			double fv = f.apply(v);
			if(fv < min) {
				min = fv;
				ret.set("x", x);
			}
			x += delta;
		}
		
		//Find out the min value on left boundary
		min = Double.MAX_VALUE;
		double y = 0.0;
		v.set("x",1.0);
		while(y < 4.0) {
			v.set("y", y);
			double fv = f.apply(v);
			if(fv < min) {
				min = fv;
				ret.set("y", y);
			}
			y += delta;
		}
		return ret;
	}

	public static void main(String[] args) {
		Curvature simulation = new Curvature();
		if(args.length == 0) {
			System.out.println("Args: <output_folder> <mesh_file> <x0> <y0> <radius> <mua> [<number of inclusions>] [<distance>]");
			System.out.println("Note:");
			System.out.println("    *distance: the distance between tow inclusions, 0.0 means just one inclusion");
			System.out.println("    *(x0,y0): the center of (first) inclusion");
			System.out.println("    *radius: the radius of each inclusion");
			System.out.println("------------Use Default Arguments-------------");
			System.out.println("output_folder=Curvature");
			System.out.println("mesh_file=prostate_test9_ex.grd");
			System.out.println("x0=2.0");
			System.out.println("y0=2.2");
			System.out.println("radius=0.5");
			System.out.println("mua=0.2");
			System.out.println("number of inclusions=1");
			System.out.println("distance=0.0");
			System.out.println("----------------------------------------------");
		}
		else{
			simulation.outputFolder = args[0];
			simulation.mesh_file = args[1];
			simulation.x0 = Double.parseDouble(args[2]);
			simulation.y0 = Double.parseDouble(args[3]);
			simulation.radius = Double.parseDouble(args[4]);
			simulation.mua = Double.parseDouble(args[5]);
			if(args.length >= 7)
				simulation.num_inclusions = Integer.parseInt(args[6]);
			if(args.length >= 8)
				simulation.distance = Double.parseDouble(args[7]);
			System.out.println("------------Use User Arguments-------------");
			System.out.println("output_folder="+simulation.outputFolder);
			System.out.println("mesh_file="+simulation.mesh_file);
			System.out.println("x0="+simulation.x0);
			System.out.println("y0="+simulation.y0);
			System.out.println("radius="+simulation.radius);
			System.out.println("mua="+simulation.mua);
			if(args.length >= 7)
				System.out.println("number of inclusions="+simulation.num_inclusions);
			else
				System.out.println("number of inclusions=1");
			if(args.length >= 8)
				System.out.println("distance="+simulation.distance);
			else
				System.out.println("distance=0.0");
			System.out.println("----------------------------------------------");
		}
		simulation.run();
	}

}
