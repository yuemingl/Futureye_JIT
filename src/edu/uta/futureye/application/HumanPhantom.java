package edu.uta.futureye.application;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.uta.futureye.algebra.SparseVectorHashMap;
import edu.uta.futureye.algebra.intf.SparseVector;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.geometry.Point;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.Vector2Function;
import edu.uta.futureye.function.operator.FMath;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.Utils;
import edu.uta.futureye.util.container.NodeList;

/**
 * 处理Fenghua数据：
 *   Research/matlab/HumanPhantom(Fenghua)
 * 
 * @author liuyueming
 *
 */
public class HumanPhantom {
	String outputFolder = "HumanPhantom";
	
	String gridFile = "human_phantom.grd";
	int NN = 20; //20*20 grid
	Mesh mesh = null;
	NodeList topList = null;
	
	ModelDOT model = new ModelDOT();

	
	/**
	 * [-2.8,2.8]*[-2.8,2.8]
	 * @param mesh
	 * @param u
	 * @return
	 */
	public Vector computeTailLeftLight(Mesh mesh, Vector u) {
		NodeList nodes = mesh.getNodeList();
		double top = 2.8;
		double left = -2.8;
		double width = 5.6;
		Point lPos = model.getLightPosition();
		Vector2Function fu = new Vector2Function(u,mesh,"x","y");
		SparseVector rlt = new SparseVectorHashMap(u.getDim());
		
		for(int i=1;i<=nodes.size();i++) {
			Node node = nodes.at(i);
			double len = Utils.computeLength(lPos, node);
			if(len>1 && len < width) {
				Variable var = new Variable("x",left+len).set("y", top);
				rlt.set(i, fu.apply(var));
			} else {
				rlt.set(i, u.get(i));
			}
		}
		return rlt;
	}
	
	/**
	 * [-2.8,2.8]*[-2.8,2.8]
	 * @param mesh
	 * @param u
	 * @return
	 */
	public Vector computeTailRightLight(Mesh mesh, Vector u) {
		NodeList nodes = mesh.getNodeList();
		double top = 2.8;
		double right = 2.8;
		double width = 5.6;
		Point lPos = model.getLightPosition();
		Vector2Function fu = new Vector2Function(u,mesh,"x","y");
		SparseVector rlt = new SparseVectorHashMap(u.getDim());
		
		for(int i=1;i<=nodes.size();i++) {
			Node node = nodes.at(i);
			double len = Utils.computeLength(lPos, node);
			if(len>1 && len < width) {
				Variable var = new Variable("x",right-len).set("y", top);
				rlt.set(i, fu.apply(var));
			} else {
				rlt.set(i, u.get(i));
			}
		}
		return rlt;
	}
	
//	double []right = {-0.00537564035466741,	
//	-0.00537564035466741,
//	-0.00537564035466741,
//	-0.00537564035466741,
//	-0.00537564035466741,
//	-0.00510100530947088,
//	-0.00400246512868475,
//	-0.00290392494789861,
//	-0.00180538476711248,
//	-0.000706844586326348,
//	-0.000322799079951417,
//	-0.000653248247987694,
//	-0.000983697416023970,
//	-0.00131414658406025,
//	-0.00164459575209652,
//	-0.00172720804410559,
//	-0.00172720804410559,
//	-0.00172720804410559,
//	-0.00172720804410559,
//	-0.00172720804410559};
	
//	double[] left ={-0.00197755193246076,
//	-0.00197755193246076,
//	-0.00197755193246076,
//	-0.00197755193246076,
//	-0.00197755193246076,
//	-0.00210946948307668,
//	-0.00263713968554036,
//	-0.00316480988800404,
//	-0.00369248009046772,
//	-0.00422015029293140,
//	-0.00545555122437696,
//	-0.00739868288480441,
//	-0.00934181454523186,
//	-0.0112849462056593,
//	-0.0132280778660867,
//	-0.0137138607811936,
//	-0.0137138607811936,
//	-0.0137138607811936,
//	-0.0137138607811936,
//	-0.0137138607811936};


	/**
	 * 
	 * @param mesh
	 * @param y
	 * @param sortAsc: sort by x
	 * @return
	 */
	public NodeList getNodeListOnY(Mesh mesh, double y,final boolean sortAsc) {
		NodeList nodes = mesh.getNodeList();
		NodeList rlt = new NodeList();
		for(int i=1;i<=nodes.size();i++) {
			Node node = nodes.at(i);
			if(Math.abs(node.coord(2)-y)<Constant.meshEps) {
				rlt.add(node);
			}
		}
		List<Node> l = rlt.toList();
		Collections.sort(l, new Comparator<Node>() {
			@Override
			public int compare(Node o1, Node o2) {
				if(o1.coord(1)>o2.coord(1))
					return sortAsc?1:-1;
				else
					return sortAsc?-1:1;
			}
		});
		return rlt;
	}
	
	public void init() {
		model.setMu_a(FC.c(0.1));
		
		MeshReader reader = new MeshReader(gridFile);
		mesh = reader.read2DMesh();		
		//Use element library to assign degree of freedom (DOF) to element
		Tools.assignLinearShapFunction(mesh);
		mesh.computeNodeBelongsToElements();
		mesh.computeNeighborNodes();
		
		topList = this.getNodeListOnY(mesh, 2.8, true);
		//System.out.println(topList);
		
	}

	public Vector run(double[] leftData, double[] rightData, int sliceNo, String outSubFloder, boolean debug) {
		String outputFolder = String.format("%s/%s/%02d",this.outputFolder,outSubFloder,sliceNo);
		String outputFolder2 = String.format("%s/%s",this.outputFolder,outSubFloder);

		int pad = 4;
	
		double ampFactor = 10;
		model.setLightPosition(-2.8, 2.8);
		Vector uLeft = model.solveNeumann(mesh);
		if(debug) Tools.plotVector(mesh, outputFolder, "uLeft.dat", uLeft);
		
		Vector tailLeftLight0 = computeTailLeftLight(mesh, uLeft);
		if(debug) Tools.plotVector(mesh, outputFolder, "tailLeftLight0.dat", tailLeftLight0);
		Vector aLeft0 = Tools.solveParamInverse(mesh, tailLeftLight0, 
				model.getDelta(), model.k,FC.c(0.1));
		if(debug) Tools.plotVector(mesh, outputFolder, "aLeft0.dat", aLeft0);

		Vector leftTest = new SparseVectorHashMap(uLeft.getDim());
		for(int i=pad;i<leftData.length-pad;i++) {
			Node node = topList.at(i+1);
			uLeft.set(node.globalIndex,uLeft.get(node.globalIndex)+ampFactor*leftData[i]);
			leftTest.set(node.globalIndex, leftData[i]);
		}
		if(debug) Tools.plotVector(mesh, outputFolder, "leftTest.dat", leftTest);
		Vector tailLeftLight = computeTailLeftLight(mesh, uLeft);
		if(debug) Tools.plotVector(mesh, outputFolder, "tailLeftLight.dat", tailLeftLight);
		Vector aLeft = Tools.solveParamInverse(mesh, tailLeftLight, 
				model.getDelta(), model.k,FC.c(0.1));
		Tools.plotVector(mesh, outputFolder, "aLeft.dat", aLeft);
		Tools.plotVector(mesh, outputFolder, "aLeftDiff.dat", aLeft.add(-1.0, aLeft0));
		
		aLeft = Utils.gaussSmooth(mesh, aLeft, 2, 0.5);
		aLeft = Utils.gaussSmooth(mesh, aLeft, 2, 0.5);
		aLeft = Utils.gaussSmooth(mesh, aLeft, 2, 0.5);
		Tools.plotVector(mesh, outputFolder, "aLeftSmooth.dat", aLeft);

		
		
		ampFactor = 10;
		model.setLightPosition(2.8, 2.8);
		Vector uRight = model.solveNeumann(mesh);
		if(debug) Tools.plotVector(mesh, outputFolder, "uRight.dat", uRight);
		
		Vector tailRightLight0 = computeTailRightLight(mesh, uRight);
		if(debug) Tools.plotVector(mesh, outputFolder, "tailRightLight0.dat", tailRightLight0);
		Vector aRight0 = Tools.solveParamInverse(mesh, tailRightLight0, 
				model.getDelta(), model.k,FC.c(0.1));
		if(debug) Tools.plotVector(mesh, outputFolder, "aRight0.dat", aRight0);
		
		Vector rightTest = new SparseVectorHashMap(uLeft.getDim());
		for(int i=pad;i<rightData.length-pad;i++) {
			Node node = topList.at(i+1);
			uRight.set(node.globalIndex,uRight.get(node.globalIndex)+ampFactor*rightData[i]);
			rightTest.set(node.globalIndex, rightData[i]);
		}
		if(debug) Tools.plotVector(mesh, outputFolder, "rightTest.dat", rightTest);
		Vector tailRightLight = computeTailRightLight(mesh, uRight);
		if(debug) Tools.plotVector(mesh, outputFolder, "tailRightLight.dat", tailRightLight);
		Vector aRight = Tools.solveParamInverse(mesh, tailRightLight, 
				model.getDelta(), model.k,FC.c(0.1));
		Tools.plotVector(mesh, outputFolder, "aRight.dat", aRight);		
		Tools.plotVector(mesh, outputFolder, "aRightDiff.dat", aRight.add(-1.0, aRight0));
		aRight = Utils.gaussSmooth(mesh, aRight, 2, 0.5);
		aRight = Utils.gaussSmooth(mesh, aRight, 2, 0.5);
		aRight = Utils.gaussSmooth(mesh, aRight, 2, 0.5);
		Tools.plotVector(mesh, outputFolder, "aRightSmooth.dat", aRight);

		Vector aAvg = FMath.axpy(1.0, aLeft, aRight);
		Tools.plotVector(mesh, outputFolder2, String.format("aAvgSmooth_%02d.dat",sliceNo), aAvg);
		
		aAvg = Utils.gaussSmooth(mesh, aAvg, 2, 0.5);
		aAvg = Utils.gaussSmooth(mesh, aAvg, 2, 0.5);
		aAvg = Utils.gaussSmooth(mesh, aAvg, 2, 0.5);
		for(int i=1;i<aAvg.getDim();i++) {
			if(aAvg.get(i)<0.0) aAvg.set(i,0.0);
		}
		Tools.plotVector(mesh, outputFolder2, String.format("aAvgSmoothCut_%02d.dat",sliceNo), aAvg);
			
		return aAvg;

	}
	
	public Mesh read3DMesh(String file) {
		MeshReader reader = new MeshReader(file);
		Mesh mesh = reader.read3DMesh();
		//Vector v = new SparseVectorHashMap(mesh.getNodeList().size());
		//Tools.plotVector(mesh, outputFolder, "3D.dat", v);
		return mesh;
	}
	
	
	/**
	 * L->left data
	 * R->right data
	 * T->top data
	 * B->bottom data
	 * @param dataFolder
	 * @return
	 */
	public Map<String,List<double[]>> readAllData(String dataFolder) {
		String[] all = {"L","R","T","B"};
		Map<String,List<double[]>> rlt = new HashMap<String,List<double[]>>();
		for(int j=0;j<all.length;j++) {
			List<double[]> list = new ArrayList<double[]>();
			for(int i=1;i<=NN;i++) {
				String file = String.format("./HumanPhantom/%s/HumanPhantom_%s%d.txt", 
						dataFolder,	all[j], i);
				double[] data = readData(file);
				list.add(data);
			}
			rlt.put(all[j],list);
		}
		return rlt;
	}

	public double [] readData(String fileName) {
		FileInputStream in;
		try {
			in = new FileInputStream(fileName);

			InputStreamReader reader = new InputStreamReader(in,"UTF-8");
			BufferedReader br = new BufferedReader(reader);
	
			String str = null;
			double[] rlt = null;
			if((str = br.readLine()) != null){
				String[] line = str.split("(\\s)+");
				rlt = new double[line.length];
				for(int i=0;i<line.length;i++)
					rlt[i] = Double.parseDouble(line[i]);
			}
			br.close();
			in.close();
			return rlt;
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}
	
	/**
	 * [-2.8,2.8]*[-2.8,2.8]*[-2.8,2.8]
	 * @param alphaList
	 */
	public void buildResult3D(List<Vector> alphaList,String rltFile) {
		Mesh mesh3D = read3DMesh("human_phantom3D.grd");
		Vector v3D = new SparseVectorHashMap(mesh3D.getNodeList().size());
		NodeList nodes = mesh3D.getNodeList();
		double dz = 5.8/NN;
		List<Vector2Function> faList = new ArrayList<Vector2Function>();
		for(int i=0;i<alphaList.size();i++) {
			faList.add(new Vector2Function(alphaList.get(i),mesh,"x","y"));
		}
		for(int i=1;i<=nodes.size();i++) {
			Node node = nodes.at(i);
			int index = (int)((node.coord(3) - (-2.8) - 0.001)/dz);
			v3D.set(i, faList.get(index).apply(
					new Variable("x",node.coord(1)).set("y",node.coord(2)))
					);
		}
		Tools.plotVector(mesh3D, outputFolder, rltFile, v3D);
	}
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		HumanPhantom hp = new HumanPhantom();
		hp.init();

		hp.plot3DInclusion();
		
		//String folder = "HumanPhantom_d30";
		//String folder = "HumanPhantom_d30d30";
		String folder = "HumanPhantom_d20d30";
		Map<String,List<double[]>> allData = hp.readAllData(folder);
		
		List<double[]> LD = allData.get("L");
		List<double[]> RD = allData.get("R");
		List<Vector> aList = new ArrayList<Vector>();
		for(int i=0;i<LD.size();i++)
			aList.add(hp.run(LD.get(i),RD.get(i),i+1,"out"+folder,false));
		
		hp.buildResult3D(aList,"3D"+folder+".dat");
	}
	
	public void plot3DInclusion() {
		Mesh mesh3D = read3DMesh("./HumanPhantom/inclusion.grd");
		Vector v3D = new SparseVectorHashMap(mesh3D.getNodeList().size());
		Tools.plotVector(mesh3D, outputFolder, "inclusion.dat", v3D);
	}

}
