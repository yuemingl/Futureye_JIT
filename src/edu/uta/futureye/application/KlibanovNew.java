package edu.uta.futureye.application;

import cern.jet.math.Constants;
import edu.uta.futureye.algebra.SpaceVector;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.function.AbstractMathFun;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FDelta;
import edu.uta.futureye.function.basic.Vector2Function;
import edu.uta.futureye.function.intf.MathFun;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.container.NodeList;
import static edu.uta.futureye.function.operator.FMath.*;

/**
 * New Algorithm for DOT in CW
 * 
 * @author yuemingl
 *
 */
public class KlibanovNew {
	String outputFolder = "KlibanovNew";
	Mesh omega = null;
	Mesh omega0 = null;
	class DomainInfo {
		double x_min = -3.5;
		double x_max = 3.5;
		double z_min = -0.5;
		double z_max = 3.5;
		
		int light_num = 10;
		double[] light_x = range(-3.5, 3.5, light_num);
		double[] light_y = range(4.0,4.0,light_num);
		
		double z1= 4.0;
	}
	DomainInfo domain = new DomainInfo();
	
	//Parameters of forward problem
	double k = 0.1;
	MathFun ax = ModelParam.getMu_a(-0.5, 2.5, 0.5, //center: (x,y;r)
									0.3, //maxMu_a: 0.3
									1); //type: one inclusion
	
	/**
	 * Load mesh
	 */
	public void init() {
		//String gridFileBig = "KlibanovNewOmega0.grd";
		//String gridFileSmall = "KlibanovNewOmega.grd";
		String gridFileBig = "backscatter1_omega0.grd";
		String gridFileSmall = "backscatter1_omega.grd";
		MeshReader readerForward = new MeshReader(gridFileBig);
		omega0 = readerForward.read2DMesh();
		MeshReader readerGCM = new MeshReader(gridFileSmall);
		omega = readerGCM.read2DMesh();
		
		//Use element library to assign degree of freedom (DOF) to element
		Tools.assignLinearShapFunction(omega0);
		Tools.assignLinearShapFunction(omega);
		omega0.computeNodeBelongsToElements();
		omega0.computeNeighborNodes();
		omega.computeNodeBelongsToElements();
		omega.computeNeighborNodes();

	}
	
	/**
	 * Laplace transform (Fourier transform)
	 * 
	 * U(x,z,s) = ReU + i*ImU = int_{-\inf}^{\inf}{u(x,z,x0)*e^{isx0}}dx0
	 * 
	 * where i=sqrt(-1)
	 * 
	 * @param u
	 * @param ReU
	 * @param ImU
	 */
	public void transformLaplace(
			Vector[] u, //input: Light Intensity
			double s,   //input: Pseudo Frequency
			Vector ReU, //output: Real Part of U
			Vector ImU  //output: Image Part of U
			) {
		ReU.setAll(0.0);
		ImU.setAll(0.0);
		int dim = u[0].getDim();
		double dx0 = domain.light_x[1] - domain.light_x[0]; 
		for(int i=0; i<domain.light_num; ++i) {
			double x0 = domain.light_x[i];
			for(int j=1; j<=dim; ++j) {
				ReU.add(j, u[i].get(j)*Math.cos(s*x0)*dx0);
				ImU.add(j, u[i].get(j)*Math.sin(s*x0)*dx0);
			}
		}
	}
	
	/**
	 * Solve a group of forward problems with different light sources alone a line 
	 * L = {x \in R, z=z1}
	 * 
	 * @return
	 */
	public Vector[] solveForwardProblem() {
		ModelDOT model = new ModelDOT();
		model.setMu_sp(10.0);
		model.setMu_a(ax.A(k)); //k+a(x) 
		Tools.plotFunction(omega0, outputFolder, "mu_a_real.dat", model.getMu_a());
		
		Vector[] ret = new Vector[domain.light_num];
		for(int i=0; i<domain.light_num; ++i) {
			model.setLightPosition(domain.light_x[i], domain.light_y[i]);
			ret[i] = model.solveNeumann(omega0);
			//Tools.plotVector(omega0, outputFolder, String.format("u%d.dat", i), ret[i]);
		}
		
		return ret;
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
	
	/**
	 * w = U*e^(-isx)
	 * @param ReU
	 * @param ImU
	 * @param Rew
	 * @param Imw
	 */
	public void compute_w(
			double s, //pseudo frequency
			Vector ReU, Vector ImU, //input
			Vector Rew, Vector Imw  //output
			) {
		NodeList nodes = omega.getNodeList();
		for(int i=1; i<=nodes.size(); i++) {
			Node n = nodes.at(i);
			Rew.set(i, ReU.get(i)*Math.cos(-s*n.coord(1)));
			Imw.set(i,ImU.get(i)*Math.sin(-s*n.coord(1)));
		}
	}
	
	/**
	 * Solve 
	 * \Delta{w} + 2isw_x - (s^2+k^+a(x))*w = -delta(z-z1)
	 * 
	 * @param s
	 * @param Rephi Boundary condition
	 * @param Imphi Boundary condition
	 */
	public void solver_w(double s, Vector Rephi, Vector Imphi) {
		ModelPoissonEx model_Rew = new ModelPoissonEx(omega);
		ModelPoissonEx model_Imw = new ModelPoissonEx(omega);
		
		//delta(z-z1) == delta(y-y1) [let y denote z]
		MathFun delta = new AbstractMathFun("x","y") {
			FDelta dt = new FDelta(new Variable("x",domain.z1),0.01,2e5);
			Variable tmp = new Variable();
			public double apply(Variable v) {
				tmp.set("x", v.get("y"));
				return dt.apply(tmp);
			}
		};
		Tools.plotFunction(omega, outputFolder, String.format("delta_z.dat"), delta);

		//Return a positive value to indicate the Dirichlet boundary type
		MathFun diriMark = new AbstractMathFun("x","y") {
			public double apply(Variable v) {
				if(Math.abs(domain.z_max-v.get("y"))<Constant.meshEps)
					return 1;
				return -1;
			}
		};
		
		model_Rew.f = delta;
		model_Rew.k = C1;
		model_Rew.c = ax.A(s*s+k*k);
		Vector Rew = model_Rew.solveMixedBorder(diriMark , new Vector2Function(Rephi,omega,"x","y"), C0, C(k));
		Tools.plotVector(omega, outputFolder, String.format("Rew_iter0.dat"), Rew);
		
		//Solve iteratively
		Vector Imw = null;
		for(int i=1; i<=10; ++i) {
			Vector Rew_x = Tools.computeDerivative(omega, Rew, "x");
			model_Imw.f = new Vector2Function(Rew_x,omega,"x","y").M(2*s);
			model_Imw.k = C1;
			model_Imw.c = ax.A(s*s+k*k);
			Imw = model_Rew.solveMixedBorder(diriMark , new Vector2Function(Imphi,omega,"x","y"), C0, C(k));
			Tools.plotVector(omega, outputFolder, String.format("Imw_iter%d.dat",i), Imw);
		
			Vector Imw_x = Tools.computeDerivative(omega, Imw, "x");
			model_Rew.f = delta.S(new Vector2Function(Imw_x,omega,"x","y").M(2*s));
			Rew = model_Rew.solveMixedBorder(diriMark , new Vector2Function(Rephi,omega,"x","y"), C0, C(k));
			Tools.plotVector(omega, outputFolder, String.format("Rew_iter%d.dat",i), Rew);
		}
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		KlibanovNew prm = new KlibanovNew();
		prm.init();
		Vector[] u0 = prm.solveForwardProblem(); //forward solution on mesh omega0
		Vector[] u = new Vector[u0.length];
		for(int i=0; i<u0.length; ++i) {
			u[i] = Tools.extractData(prm.omega0, prm.omega, u0[i]); //forward solution on mesh omega
		}
		
		int dim = u[0].getDim();
		Vector ReU = new SpaceVector(dim);
		Vector ImU = new SpaceVector(dim);
		double s = 5.0;// pseudo frequency
		prm.transformLaplace(u,s,ReU,ImU);
		Tools.plotVector(prm.omega, prm.outputFolder, String.format("ReU.dat"), ReU);
		Tools.plotVector(prm.omega, prm.outputFolder, String.format("ImU.dat"), ImU);
		
		Vector Rew = new SpaceVector(dim);
		Vector Imw = new SpaceVector(dim);
		prm.compute_w(s, ReU, ImU, Rew, Imw);
		Tools.plotVector(prm.omega, prm.outputFolder, String.format("Rew.dat"), Rew);
		Tools.plotVector(prm.omega, prm.outputFolder, String.format("Imw.dat"), Imw);
		
		Vector Rew_x = Tools.computeDerivative(prm.omega, Rew, "x");
		Vector Imw_x = Tools.computeDerivative(prm.omega, Imw, "x");
		Tools.plotVector(prm.omega, prm.outputFolder, String.format("Rew_x.dat"), Rew_x);
		Tools.plotVector(prm.omega, prm.outputFolder, String.format("Imw_x.dat"), Imw_x);
		
		//Solve w iteratively
		prm.solver_w(s, Rew, Imw);
		
		//update tail
		//write you code here to update tail
		
		//solve q
		//Please modify the file GCMModelNew.java to implement 
		//the new algorithm based on the current implementation of GCM
		//GCMModelNew model = new GCMModelNew(prm.outputFolder);
		//model.solveGCM(mesh, N, s, phi, tailT)
	}

}
