package edu.uta.futureye.application;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;

import edu.uta.futureye.algebra.SpaceVector;
import edu.uta.futureye.algebra.SparseVectorHashMap;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.algebra.solver.external.SolverJBLAS;
import edu.uta.futureye.application.HumanReal.Part;
import edu.uta.futureye.core.EdgeLocal;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.FaceLocal;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.NodeLocal;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.core.geometry.GeoEntity3D;
import edu.uta.futureye.core.intf.Assembler;
import edu.uta.futureye.function.MultiVarFunc;
import edu.uta.futureye.function.FMath;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.FDelta;
import edu.uta.futureye.function.basic.Vector2MathFunc;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.io.MeshWriter;
import edu.uta.futureye.lib.assembler.AssemblerScalar;
import edu.uta.futureye.lib.assembler.AssemblerScalarFast;
import edu.uta.futureye.lib.element.FETrilinearHexahedron;
import edu.uta.futureye.lib.weakform.WeakFormLaplace;
import edu.uta.futureye.lib.weakform.WeakFormLaplace3D;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.NodeList;
import edu.uta.futureye.util.container.ObjList;
import static edu.uta.futureye.function.FMath.*;

public class HumanRealSimulate {
	String workFolder = "HumanReal/Simulate";
	//light source position
	public Variable lightPosition = null; 
	MathFunc delta = null;
	Mesh meshOmega = null;
	
	double mu_sp = 0.0;
	MathFunc mu_a = null;
	
	public void setModelParam(double mu_sp, MathFunc mu_a) {
		this.mu_sp = mu_sp;
		this.mu_a = mu_a;
	}
	public void setLightPosition(double x,double y,double z) {
		this.lightPosition = new Variable();
		this.lightPosition.set("x", x);
		this.lightPosition.set("y", y);
		this.lightPosition.set("z", z);
		delta = new FDelta(this.lightPosition,0.01,2e5);
	}
	public void setLightPosition(Variable v) {
		this.lightPosition = v;
		delta = new FDelta(this.lightPosition,0.01,2e5);
	}
	
	public void readMesh(String meshName) {
		MeshReader reader = new MeshReader(workFolder+"/"+meshName+".grd");
		meshOmega = reader.read3DMesh(); //3D
		meshOmega.setName(meshName);
		meshOmega.computeNodeBelongsToElements(); //worked in 3D
		meshOmega.computeGlobalEdge();
        ElementList eList = meshOmega.getElementList();
        FETrilinearHexahedron feTLH = new FETrilinearHexahedron();
        for(int i=1;i<=eList.size();i++)
        	feTLH.assignTo(eList.at(i));		
	}
	
	/**
	 * Solver the following model problem:
	 * 
	 *   -\nabla{(1/k)*\nabla{u}} + a*u = f
	 * 
	 * where
	 *   f = delta(x-x0)*delta(y-y0)*delta(z-z0)
	 *   u = u(x,y,z)
	 *   k = 1/(3*mu_s')
	 *   a = a(x,y,z) = mu_a(x,y,z)
	 */
	public Vector solveForward(String info) {
		System.out.println("Solve forward: mesh \""+meshOmega.getName() + "\" info:" + info);
				
		HashMap<NodeType, MathFunc> mapNTF = new HashMap<NodeType, MathFunc>();
		mapNTF.put(NodeType.Neumann, null);
		meshOmega.markBorderNode(mapNTF);
		meshOmega.writeNodesInfo(workFolder+"/meshInfo.dat");

//		int count=0,count2=0;
//		for(Element e : mesh.getElementList()) {
//			if(e.isBorderElement()) count++;
//			GeoEntity3D entity = e.getGeoEntity3D();
//				ObjList<FaceLocal> faces = entity.getFaces();
//				for(int i=1;i<=faces.size();i++) {
//					if(faces.at(i).isBorderFace())
//						count2++;
//				}
//			//检查边界单元的结点是否在同一个面上：有一个坐标值全相同（因为网格是正规的）
//			ElementList beList = e.getBorderElements();
//			for(Element be : beList) {
//				boolean find2 = false;
//				for(int j=1;j<=3;j++) {
//					boolean find = true;
//					for(int i=2;i<=be.nodes.size();i++) {
//						if(be.nodes.at(1).coord(j)!=be.nodes.at(i).coord(j)) {
//							find = false;
//							break;
//						}
//					}
//					if(find) { find2 = true; break; }
//				}
//				if(!find2) {
//					System.out.println("Error : "+be);
//				}
//			}
//		}
//		System.out.println(count);
//		System.out.println(count2);
		
		//User defined weak form of PDE (including bounder conditions)
		WeakFormLaplace3D weakForm = new WeakFormLaplace3D();
		
		weakForm.setF(delta);
		//weakForm.setF(C(1));
		
		weakForm.setParam(
				C(1./(3*this.mu_sp)), 
				this.mu_a, 
				null, 
				C(1./(3*this.mu_sp))//null
			);
		
		Assembler assembler = new AssemblerScalar(meshOmega, weakForm);
		System.out.println("Begin Assemble...");
		long begin = System.currentTimeMillis();
		assembler.assemble();
		Matrix stiff = assembler.getStiffnessMatrix();
		Vector load = assembler.getLoadVector();
		System.out.println("Assemble done!");
		long end = System.currentTimeMillis();
		System.out.println("Time used:"+(end-begin));
		System.out.println("imposeDirichletCondition()...");
		begin = System.currentTimeMillis();
		assembler.imposeDirichletCondition(FC.C0);
		end = System.currentTimeMillis();
		System.out.println("imposeDirichletCondition() time used:"+(end-begin));
		
		SolverJBLAS solver = new SolverJBLAS();
		begin = System.currentTimeMillis();
		Vector u = solver.solveDGESV(stiff, load);
		end = System.currentTimeMillis();
		System.out.println("SolverJBLAS time used:"+(end-begin));
	    
	    //Tools.plotVector(meshOmega, workFolder, meshName+"_u.dat", u);
		return u;
	}
	
	/**
	 * Measure surface is at y=0 of Omega
	 * Measure surface(xi,yi) <==> Omega(yi,0,xi)
	 * 
	 * @param measureSurface
	 * @param u
	 * @return
	 */
	public Vector extractMeasureSurfaceValues(Mesh measureSurface, Vector u) {
		NodeList surfNodes = measureSurface.getNodeList();
		Vector rlt = new SpaceVector(surfNodes.size());
		double []coords = new double[3];
		for(int i=1;i<=surfNodes.size();i++) {
			Node sn = surfNodes.at(i);
			coords[0] = sn.coord(2);
			coords[1] = 0.0;
			coords[2] = sn.coord(1);
			Node on = this.meshOmega.findNode(new Node().set(0,coords));
			if(on != null) 
				rlt.set(i, u.get(on.globalIndex));
		}
		return rlt;
	}
	
	public void simulateForward(
			Variable[] lightSources,
			Vector u0, //output
			Vector u //output
			) {
		//mu_s'
		double mu_sp = 10;
		//mu_a
		MathFunc fmu_a = new MultiVarFunc("x","y","z"){
			double x0 = 4.0,y0=-3.0, z0=4.0;//z0=13
			@Override
			public double apply(Variable v) {
				double x = v.get("x");
				double y = v.get("y");
				double z = v.get("z");
				if(Math.sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0))<1.0)
					return 0.5;
				else
					return 0.1;
			}
		};
		Tools.plotFunction(meshOmega, workFolder, "mu_a.dat", fmu_a);
		MathFunc fmu_a0 = C(0.1);

		
		for(int i=0;i<lightSources.length;i++) {
			setLightPosition(lightSources[i]);
			Tools.plotFunction(meshOmega, workFolder, String.format("delta_src%02d.dat",i), this.delta);
			
			setModelParam(mu_sp, fmu_a0);
			Vector _u0 = solveForward("background");
			if(i==0) u0.set(_u0); else u0.axpy(1.0, _u0);
			Tools.plotVector(meshOmega, workFolder, String.format(meshOmega.getName()+"_u0_src%02d.dat",i), _u0);
	
			setModelParam(mu_sp, fmu_a);
			Vector _u = solveForward("inclusion");
			if(i==0) u.set(_u); else u.axpy(1.0, _u);
			Tools.plotVector(meshOmega, workFolder, String.format(meshOmega.getName()+"_u_src%02d.dat",i), _u);
		}
		u0.scale(1.0/lightSources.length);
		u.scale(1.0/lightSources.length);
	}
	
	/**
	 * 将模拟Detector数据写入文件detector.txt
	 * 
	 * @param meshMeasureSurface
	 * @param u0
	 * @param u
	 * @param meaSurf0
	 * @param meaSurf
	 */
	public void simulateDetector(
			Mesh meshMeasureSurface, Vector u0, Vector u, //input
			Vector meaSurf0, Vector meaSurf //output
			) {
		Vector _meaSurf0 = extractMeasureSurfaceValues(meshMeasureSurface, u0);
		meaSurf0.set(_meaSurf0);
		Tools.plotVector(meshMeasureSurface, workFolder, meshOmega.getName()+"_u0_MeaSurf.dat", meaSurf0);
		
		Vector _meaSurf = extractMeasureSurfaceValues(meshMeasureSurface, u);
		meaSurf.set(_meaSurf);
		Tools.plotVector(meshMeasureSurface, workFolder, meshOmega.getName()+"_u_MeaSurf.dat", meaSurf);
		
		//Coordinates of detectors
		double []setup_x = {1,3,5,7,  2,  4,  6,1,3,5,7};
		double []setup_y = {3,3,3,3,4.5,4.5,4.5,6,6,6,6};
		Vector2MathFunc fmeaSurf0 = new Vector2MathFunc(meaSurf0, meshMeasureSurface, "x","y");
		Vector2MathFunc fmeaSurf = new Vector2MathFunc(meaSurf, meshMeasureSurface, "x","y");
		Vector diffIntensity = new SpaceVector(setup_x.length);
		Vector diffOD = new SpaceVector(setup_x.length);
		for(int i=0;i<setup_x.length;i++) {
			Variable v = new Variable("x",setup_x[i]).set("y",setup_y[i]);
			//diffIntensity = I - I0
			diffIntensity.set(i+1, fmeaSurf.apply(v)-fmeaSurf0.apply(v));
			//diffOD = - log (I/I0)
			diffOD.set(i+1, -Math.log(fmeaSurf0.apply(v)/fmeaSurf.apply(v)));
		}
		
		MeshWriter.writeTechplotLine(workFolder + "/detector.txt", 
				new SpaceVector(setup_x), 
				new SpaceVector(setup_y),
				diffIntensity,
				diffOD
				);	
	}
	
	/**
	 * 测试inclusion深度变化和大小变化对测量数据的影响。
	 * 主要测试一个大一些的inclusion在较深位置与一个小一些的inclusion在较浅位置
	 * 在测量面得到的光强变化如果在最暗处相等，是否能够区分这两种情况？
	 * 
	 */
	public static void testDepthBig() {
		HumanRealSimulate sim = new HumanRealSimulate();
		sim.workFolder = "HumanReal/SimulateDepth";
		String meshName = "mesh3DLeft";
		//String meshName = "mesh3DRight";
		sim.readMesh(meshName);
		
		String meshMeasureSurfaceName = "meshMeasureSurfaceLeft";
		//String meshMeasureSurfaceName = "meshMeasureSurfaceRight";
		MeshReader reader = new MeshReader(sim.workFolder+"/"+meshMeasureSurfaceName+".grd");
		Mesh meshMeasureSurface = reader.read2DMesh(); //2D
		meshMeasureSurface.setName(meshMeasureSurfaceName);
		
		Variable[] lightSources = new Variable[3];
		lightSources[0] = new Variable("x",1.5).set("y", -0.5).set("z",4);
		
		//mu_s'
		double mu_sp = 10;
		MathFunc fmu_a0 = C(0.1);
		sim.setLightPosition(lightSources[0]);
		Tools.plotFunction(sim.meshOmega, sim.workFolder, String.format("delta_src%02d.dat",0), sim.delta);
		sim.setModelParam(mu_sp, fmu_a0);
		Vector u0 = sim.solveForward("simulation background");
		Tools.plotVector(sim.meshOmega, sim.workFolder, String.format(sim.meshOmega.getName()+"_u0_src%02d.dat",0), u0);
		Vector meaSurf0 = sim.extractMeasureSurfaceValues(meshMeasureSurface, u0);
		Tools.plotVector(meshMeasureSurface, sim.workFolder, sim.meshOmega.getName()+"_u0_MeaSurf.dat", meaSurf0);

		
		//mu_a
		MathFunc fmu_aBig = new MultiVarFunc("x","y","z"){
			//深度 y=-3.0
			double x0 = 4.0,y0=-3.0, z0=4.0;//z0=13
			@Override
			public double apply(Variable v) {
				double x = v.get("x");
				double y = v.get("y");
				double z = v.get("z");
				if(Math.sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0))<1.0)
					return 0.5;
				else
					return 0.1;
			}
		};
		Tools.plotFunction(sim.meshOmega, sim.workFolder, "mu_a_big.dat", fmu_aBig);
		sim.setModelParam(mu_sp, fmu_aBig);
		Vector uBig = sim.solveForward("simulation big inclusion");
		Tools.plotVector(sim.meshOmega, sim.workFolder, String.format(sim.meshOmega.getName()+"_u_big_src%02d.dat",0), uBig);
		Vector meaSurfBig = sim.extractMeasureSurfaceValues(meshMeasureSurface, uBig);
		Tools.plotVector(meshMeasureSurface, sim.workFolder, sim.meshOmega.getName()+"_u_big_MeaSurf.dat", meaSurfBig);
		
	}
	
	public static void testDepthSmall() {
		HumanRealSimulate sim = new HumanRealSimulate();
		sim.workFolder = "HumanReal/SimulateDepth";
		String meshName = "mesh3DLeft";
		//String meshName = "mesh3DRight";
		sim.readMesh(meshName);
		
		String meshMeasureSurfaceName = "meshMeasureSurfaceLeft";
		//String meshMeasureSurfaceName = "meshMeasureSurfaceRight";
		MeshReader reader = new MeshReader(sim.workFolder+"/"+meshMeasureSurfaceName+".grd");
		Mesh meshMeasureSurface = reader.read2DMesh(); //2D
		meshMeasureSurface.setName(meshMeasureSurfaceName);
		
		Variable[] lightSources = new Variable[3];
		lightSources[0] = new Variable("x",1.5).set("y", -0.5).set("z",4);
		sim.setLightPosition(lightSources[0]);
	
		//mu_s'
		double mu_sp = 10;
		MathFunc fmu_aSmall = new MultiVarFunc("x","y","z"){
			//深度 y=-1.5
			double x0 = 4.0,y0=-1.5, z0=4.0;//z0=13
			@Override
			public double apply(Variable v) {
				double x = v.get("x");
				double y = v.get("y");
				double z = v.get("z");
				if(Math.sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0))<0.5)
					return 0.18;
				else
					return 0.1;
			}
		};
		Tools.plotFunction(sim.meshOmega, sim.workFolder, "mu_a_small.dat", fmu_aSmall);
		sim.setModelParam(mu_sp, fmu_aSmall);
		Vector uSmall = sim.solveForward("simulation small inclusion");
		Tools.plotVector(sim.meshOmega, sim.workFolder, String.format(sim.meshOmega.getName()+"_u_small_src%02d.dat",0), uSmall);
		Vector meaSurfSmall = sim.extractMeasureSurfaceValues(meshMeasureSurface, uSmall);
		Tools.plotVector(meshMeasureSurface, sim.workFolder, sim.meshOmega.getName()+"_u_small_MeaSurf.dat", meaSurfSmall);
		
	}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		//testDepthBig();
		testDepthSmall();

//		
//		HumanRealSimulate sim = new HumanRealSimulate();
//		sim.plot3DInclusion();
//
//		String meshName = "mesh3DLeft";
//		//String meshName = "mesh3DRight";
//		sim.readMesh(meshName);
//		
//		String meshMeasureSurfaceName = "meshMeasureSurfaceLeft";
//		//String meshMeasureSurfaceName = "meshMeasureSurfaceRight";
//		MeshReader reader = new MeshReader(sim.workFolder+"/"+meshMeasureSurfaceName+".grd");
//		Mesh meshMeasureSurface = reader.read2DMesh(); //2D
//		meshMeasureSurface.setName(meshMeasureSurfaceName);
//		
//		
//		Vector u0Sum = new SparseVectorHashMap();
//		Vector uSum = new SparseVectorHashMap();
//		Variable[] lightSources = new Variable[3];
//		lightSources[0] = new Variable("x",1.5).set("y", -0.5).set("z",2);
//		lightSources[1] = new Variable("x",1.5).set("y", -0.5).set("z",4);
//		lightSources[2] = new Variable("x",1.5).set("y", -0.5).set("z",6);
//		sim.simulateForward(lightSources, u0Sum, uSum);
//		
//		Vector meaSurf0 = new SparseVectorHashMap();
//		Vector meaSurf = new SparseVectorHashMap();
//		sim.simulateDetector(meshMeasureSurface, u0Sum, uSum, meaSurf0, meaSurf);
//		
//		//使用matlab进行插值 interplate_measurement_surface.m
//		
//		Vector interpMeaSurf = Tools.read2DFunctionValues(meshMeasureSurface,sim.workFolder+"/interplate_data.mat","XXI","YYI","VVI");
//		Tools.plotVector(meshMeasureSurface, sim.workFolder, meshName+"_u_MeaSurf_interp_dOD.dat", interpMeaSurf);
//		
//		Vector diffIntensity = axpy(
//				-1.0,
//				meaSurf0, 
//				axMuly(1.0, meaSurf0,exp(interpMeaSurf))
//			);
//		Tools.plotVector(meshMeasureSurface, sim.workFolder, meshName+"_u_MeaSurf_interp_dI.dat", diffIntensity);
//		Vector diffIntensityReal = axpy(-1.0,meaSurf0,meaSurf);
//		
//		
//		HumanReal hr = new HumanReal();
//		hr.init();
//		NodeList slice1Left = null;
//		NodeList slice1Right = null;
//		double[] leftData = new double[hr.NN];
//		double[] rightData = new double[hr.NN];
//		int timeIndex = 1;
//		boolean debug = true;
//		double[] leftLightCoord = {1.5,-0.5};
//		double[] rightLightCoord = {7.5,-0.5};
//		
//		//模拟插值数据
//		String outSubFloder = "simulate_interp";
//		List<Vector> alphaSlices = new ArrayList<Vector>();
//		for(int slice=1;slice<=hr.NSlice;slice++) {
//			double x = (slice-1)*0.5;
//			slice1Left = Tools.getNodeListOnX(meshMeasureSurface, x, false);
//			slice1Right = Tools.getNodeListOnX(meshMeasureSurface, x, true);
//			
//			for(int i=1;i<=slice1Left.size();i++) {
//				leftData[i-1] = diffIntensity.get(slice1Left.at(i).globalIndex);
//			}
//			for(int i=1;i<=slice1Right.size();i++) {
//				rightData[i-1] = diffIntensity.get(slice1Right.at(i).globalIndex);
//			}
//			Vector vSlice = hr.run( leftData, rightData, 
//					slice, timeIndex, outSubFloder, debug,
//					1,leftLightCoord,rightLightCoord);
//			alphaSlices.add(vSlice);
//		}
//		if(meshName.equals("mesh3DLeft"))
//			hr.buildResult3DLeft(alphaSlices, String.format("3D_Left_simulate_interp.dat"));
//		else
//			hr.buildResult3DRight(alphaSlices, String.format("3D_Right_simulate_interp.dat"));
//		
//		
//		//模拟真实数据
//		outSubFloder = "simulate_real";
//		alphaSlices.clear();
//		for(int slice=1;slice<=hr.NSlice;slice++) {
//			double x = (slice-1)*0.5;
//			slice1Left = Tools.getNodeListOnX(meshMeasureSurface, x, false);
//			slice1Right = Tools.getNodeListOnX(meshMeasureSurface, x, true);
//			
//			for(int i=1;i<=slice1Left.size();i++) {
//				leftData[i-1] = diffIntensityReal.get(slice1Left.at(i).globalIndex);
//			}
//			for(int i=1;i<=slice1Right.size();i++) {
//				rightData[i-1] = diffIntensityReal.get(slice1Right.at(i).globalIndex);
//			}
//			Vector vSlice = hr.run( leftData, rightData, 
//					slice, timeIndex, outSubFloder, debug,
//					1,leftLightCoord,rightLightCoord);
//			alphaSlices.add(vSlice);
//		}
//		if(meshName.equals("mesh3DLeft"))
//			hr.buildResult3DLeft(alphaSlices, String.format("3D_Left_simulate_real.dat"));
//		else
//			hr.buildResult3DRight(alphaSlices, String.format("3D_Right_simulate_real.dat"));
//		
	}
	
	public Mesh read3DMesh(String file) {
		MeshReader reader = new MeshReader(file);
		Mesh mesh = reader.read3DMesh();
		//Vector v = new SparseVectorHashMap(mesh.getNodeList().size());
		//Tools.plotVector(mesh, outputFolder, "3D.dat", v);
		return mesh;
	}
	
	public void plot3DInclusion() {
		Mesh mesh3D = read3DMesh(this.workFolder+"/inclusion.grd");
		Vector v3D = new SparseVectorHashMap(mesh3D.getNodeList().size());
		Tools.plotVector(mesh3D, this.workFolder, "inclusion.dat", v3D);
	}
}
