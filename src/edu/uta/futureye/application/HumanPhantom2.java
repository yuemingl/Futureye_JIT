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
import edu.uta.futureye.function.FMath;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.Vector2MathFunc;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.util.Constant;
import edu.uta.futureye.util.Utils;
import edu.uta.futureye.util.container.NodeList;

/**
 * 处理Fenghua数据：
 *   Research/matlab/HumanPhantom2(Fenghua)
 * 
 * @author liuyueming
 *
 */
public class HumanPhantom2 {
	String outputFolder = "HumanPhantom2";
	
	String gridFile = "human_phantom2.grd";
	int MM = 21;//16; //16*21 grid
	int NN = 21;
	double meshTop = 0.0;
	double meshBottom = 6.0;
	double meshLeft = -3.8;
	double meshRight = 3.8;
	Mesh mesh = null;
	NodeList topList = null;
	
	ModelDOT model = new ModelDOT();

	
	/**
	 * [-3.8,3.8]*[0,-6]
	 * @param mesh
	 * @param u
	 * @return
	 */
	public Vector computeTailLeftLight(Mesh mesh, Vector u) {
		NodeList nodes = mesh.getNodeList();
		Point rPos = model.getLightPosition();
		Vector2MathFunc fu = new Vector2MathFunc(u,mesh,"x","y");
		SparseVector rlt = new SparseVectorHashMap(u.getDim());
		double skipLen = 1.0;
		for(int i=1;i<=nodes.size();i++) {
			Node node = nodes.at(i);
			//计算当前点和光源的距离
			double len = Utils.computeLength(rPos, node);
			double dx = Math.signum(node.coord(1)-rPos.coord(1))*len*1.0;//!!!
			
			if(len>skipLen && rPos.coord(1)+dx>=meshLeft && rPos.coord(1)+dx<=meshRight) {
				Variable var = new Variable("x",rPos.coord(1)+dx).set("y", meshTop);
				rlt.set(i, fu.apply(var));
			} else {
				rlt.set(i, u.get(i));				
			}
		}
		return rlt;
	}
	
	/**
	 * [-3.8,3.8]*[0,-6]
	 * @param mesh
	 * @param u
	 * @return
	 */
	public Vector computeTailRightLight(Mesh mesh, Vector u) {
		return computeTailLeftLight(mesh,u);
	}
	
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
		
		MeshReader reader = new MeshReader(this.outputFolder+"/"+gridFile);
		mesh = reader.read2DMesh();		
		//Use element library to assign degree of freedom (DOF) to element
		Tools.assignLinearShapFunction(mesh);
		mesh.computeNodeBelongsToElements();
		mesh.computeNeighborNodes();
		
		topList = this.getNodeListOnY(mesh, this.meshTop, true);
		//System.out.println(topList);
		
	}
	
	/**
	 * Scale v1,v2 to the same scale
	 * @param v1
	 * @param v2
	 * @param out1
	 * @param out2
	 */
	public void scale(Vector v1, Vector v2, Vector out1, Vector out2) {
		double max1 = FMath.max(v1);
		double min1 = FMath.min(v1);
		double c1 = (max1 + min1)/2.0;
		double max2 = FMath.max(v2);
		double min2 = FMath.min(v2);
		double c2 = (max2 + min2)/2.0;
		//按照最大的scale有问题，改为按照v2来scale
		if(Math.abs(max1 - c1) > Constant.eps && Math.abs(c1-min1) > Constant.eps) {
			out1.set(v1);
			out1.shift(-c1);
			for(int i=1;i<=out1.getDim();i++) {
				if(out1.get(i)>=0) {
					out1.set(i,out1.get(i)*(max2-c2)/(max1-c1));
				} else {
					out1.set(i,out1.get(i)*(c2-min2)/(c1-min1));
				}
			}
			out1.shift(c2);
		} else {
			out1.set(v1);
		}
		out2.set(v2);
	}

	/**
	 * 先合并tail再计算反问题
	 * @param leftData
	 * @param rightData
	 * @param sliceNo
	 * @param outSubFloder
	 * @param debug
	 * @param ampFactor
	 * @param leftLightCoord
	 * @param rightLightCoord
	 * @return
	 */
	public Vector run1(double[] leftData, double[] rightData, int sliceNo, String outSubFloder, boolean debug,
			double ampFactor,double[] leftLightCoord, double[] rightLightCoord) {

		String outputFolder = String.format("%s/%s/%02d",this.outputFolder,outSubFloder,sliceNo);
		String outputFolder2 = String.format("%s/%s",this.outputFolder,outSubFloder);
		String timeFormat = String.format("_t%05d", 0);

		int padLeft = 6;
		int padRight1 = 6;
		int padRight2 = 6;
		////double ampFactor = 1000;
		model.setLightPosition(leftLightCoord[0], leftLightCoord[1]);
		Vector uLeft = model.solveNeumann(mesh);
		if(debug) Tools.plotVector(mesh, outputFolder, "uLeft"+timeFormat+".dat", uLeft);

		Vector tailLeftLight0 = computeTailLeftLight(mesh, uLeft);
		if(debug) Tools.plotVector(mesh, outputFolder, "tailLeftLight0"+timeFormat+".dat", tailLeftLight0);
		Vector aLeft0 = Tools.solveParamInverse(mesh, tailLeftLight0, 
				model.getDelta(), model.k,FC.c(0.1));
		if(debug) Tools.plotVector(mesh, outputFolder, "aLeft0"+timeFormat+".dat", aLeft0);

		Vector leftRawData = new SparseVectorHashMap(uLeft.getDim());
		Vector leftRawDataFull = new SparseVectorHashMap(uLeft.getDim());
		for(int i=0;i<leftData.length;i++) {
			Node node = topList.at(i+1);
			if(i>=padLeft && i<leftData.length-padLeft) {
				uLeft.set(node.globalIndex,uLeft.get(node.globalIndex)+ampFactor*leftData[i]);
				leftRawData.set(node.globalIndex, leftData[i]);
			}
			leftRawDataFull.set(node.globalIndex, leftData[i]);
		}
		if(debug) Tools.plotVector(mesh, outputFolder, "leftRawData"+timeFormat+".dat", leftRawData);
		if(debug) Tools.plotVector(mesh, outputFolder, "leftRawDataFull"+timeFormat+".dat", leftRawDataFull);
		
		
		//Left tail
		Vector tailLeftLight = computeTailLeftLight(mesh, uLeft);
		Vector tailLeftLightDiff = FMath.axpy(-1, tailLeftLight0,tailLeftLight);
		if(debug) {
			Tools.plotVector(mesh, outputFolder, "tailLeftLight"+timeFormat+".dat", tailLeftLight);
			Tools.plotVector(mesh, outputFolder, "tailLeftLightDiff"+timeFormat+".dat", tailLeftLightDiff);
		}

		model.setLightPosition(rightLightCoord[0], rightLightCoord[1]);
		Vector uRight = model.solveNeumann(mesh);
		if(debug) Tools.plotVector(mesh, outputFolder, "uRight"+timeFormat+".dat", uRight);
		
		Vector tailRightLight0 = computeTailRightLight(mesh, uRight);
		if(debug) Tools.plotVector(mesh, outputFolder, "tailRightLight0"+timeFormat+".dat", tailRightLight0);
		Vector aRight0 = Tools.solveParamInverse(mesh, tailRightLight0, 
				model.getDelta(), model.k,FC.c(0.1));
		if(debug) Tools.plotVector(mesh, outputFolder, "aRight0"+timeFormat+".dat", aRight0);
		
		Vector rightRawData = new SparseVectorHashMap(uRight.getDim());
		Vector rightRawDataFull = new SparseVectorHashMap(uRight.getDim());
		for(int i=0;i<rightData.length;i++) {
			Node node = topList.at(i+1);
			if(i>=padRight1 && i<rightData.length-padRight2) {
				uRight.set(node.globalIndex,uRight.get(node.globalIndex)+ampFactor*rightData[i]);
				rightRawData.set(node.globalIndex, rightData[i]);
			}
			rightRawDataFull.set(node.globalIndex, rightData[i]);
		}
		if(debug) Tools.plotVector(mesh, outputFolder, "rightRawData"+timeFormat+".dat", rightRawData);
		if(debug) Tools.plotVector(mesh, outputFolder, "rightRawDataFull"+timeFormat+".dat", rightRawDataFull);

		//Right tail
		Vector tailRightLight = computeTailRightLight(mesh, uRight);
		Vector tailRightLightDiff = FMath.axpy(-1, tailRightLight0,tailRightLight);
		if(debug) {
			Tools.plotVector(mesh, outputFolder, "tailRightLight"+timeFormat+".dat", tailRightLight);
			Tools.plotVector(mesh, outputFolder, "tailRightLightDiff"+timeFormat+".dat", tailRightLightDiff);
		}

		Vector tailRightLightDiffScale = tailRightLightDiff.copy();
		Vector tailLeftLightDiffScale = tailRightLightDiff.copy();
		scale(tailLeftLightDiff,tailRightLightDiff,tailLeftLightDiffScale,tailRightLightDiffScale);
		if(debug) {
			Tools.plotVector(mesh, outputFolder, "tailLeftLightDiffScale"+timeFormat+".dat", tailLeftLightDiffScale);
			Tools.plotVector(mesh, outputFolder, "tailRightLightDiffScale"+timeFormat+".dat", tailRightLightDiffScale);
		}
	
		Vector avgTailDiff = FMath.axpy(1.0, tailLeftLightDiffScale, tailRightLightDiffScale);
		//cut 70%
		double min = FMath.min(avgTailDiff);
		for(int i=1;i<avgTailDiff.getDim();i++) {
			if(avgTailDiff.get(i) > min*0.8) avgTailDiff.set(i,0.0);
		}
		if(debug) {
			Tools.plotVector(mesh, outputFolder, "tailAvgDiff_Cut"+timeFormat+".dat", avgTailDiff);
		}
	
		//只需要计算一个反问题
		Vector aRight = Tools.solveParamInverse(mesh, FMath.axpy(1.0, tailRightLight0,avgTailDiff), 
				model.getDelta(), model.k,FC.c(0.1));
		Tools.plotVector(mesh, outputFolder, "aRight"+timeFormat+".dat", aRight);		
		Tools.plotVector(mesh, outputFolder, "aRightDiff"+timeFormat+".dat", aRight.add(-1.0, aRight0));
		//去掉负值，再光滑
		aRight = Utils.gaussSmooth(mesh, aRight, 1, 0.5);
		aRight = Utils.gaussSmooth(mesh, aRight, 1, 0.4);
		aRight = Utils.gaussSmooth(mesh, aRight, 1, 0.3);
		aRight = Utils.gaussSmooth(mesh, aRight, 2, 0.4);
		aRight = Utils.gaussSmooth(mesh, aRight, 2, 0.3);
		for(int i=1;i<aRight.getDim();i++) {
			if(aRight.get(i)<0.0) aRight.set(i,0.0);
		}
		Tools.plotVector(mesh, outputFolder, "aRightSmooth_Cut"+timeFormat+".dat", aRight);

		Vector aAvg = aRight;
		
		Tools.plotVector(mesh, outputFolder2, String.format("aAvgSmooth_%02d"+timeFormat+".dat",sliceNo), aAvg);
		
		return aAvg;		
	}

	
	public Vector run(double[] leftData, double[] rightData, int sliceNo, String outSubFloder, boolean debug,
			double ampFactor,double[] leftLightCoord, double[] rightLightCoord) {

		String outputFolder = String.format("%s/%s/%02d",this.outputFolder,outSubFloder,sliceNo);
		String outputFolder2 = String.format("%s/%s",this.outputFolder,outSubFloder);
		String timeFormat = String.format("_t%05d", 0);

		int padLeft = 7;
		int padRight1 = 7;
		int padRight2 = 7;
		////double ampFactor = 1000;
		model.setLightPosition(leftLightCoord[0], leftLightCoord[1]);
		Vector uLeft = model.solveNeumann(mesh);
		if(debug) Tools.plotVector(mesh, outputFolder, "uLeft"+timeFormat+".dat", uLeft);

		Vector tailLeftLight0 = computeTailLeftLight(mesh, uLeft);
		if(debug) Tools.plotVector(mesh, outputFolder, "tailLeftLight0"+timeFormat+".dat", tailLeftLight0);
		Vector aLeft0 = Tools.solveParamInverse(mesh, tailLeftLight0, 
				model.getDelta(), model.k,FC.c(0.1));
		if(debug) Tools.plotVector(mesh, outputFolder, "aLeft0"+timeFormat+".dat", aLeft0);

		Vector leftRawData = new SparseVectorHashMap(uLeft.getDim());
		Vector leftRawDataFull = new SparseVectorHashMap(uLeft.getDim());
		for(int i=0;i<leftData.length;i++) {
			Node node = topList.at(i+1);
			if(i>=padLeft && i<leftData.length-padLeft) {
				uLeft.set(node.globalIndex,uLeft.get(node.globalIndex)+ampFactor*leftData[i]);
				leftRawData.set(node.globalIndex, leftData[i]);
			}
			leftRawDataFull.set(node.globalIndex, leftData[i]);
		}
		if(debug) Tools.plotVector(mesh, outputFolder, "leftRawData"+timeFormat+".dat", leftRawData);
		if(debug) Tools.plotVector(mesh, outputFolder, "leftRawDataFull"+timeFormat+".dat", leftRawDataFull);
		
		
		//Left tail
		Vector tailLeftLight = computeTailLeftLight(mesh, uLeft);
		if(debug) {
			Tools.plotVector(mesh, outputFolder, "tailLeftLight"+timeFormat+".dat", tailLeftLight);
			Tools.plotVector(mesh, outputFolder, "tailLeftLightDiff"+timeFormat+".dat", 
					FMath.axpy(-1, tailLeftLight0,tailLeftLight));
		}
		//smooth tail
		//tailLeftLight = Utils.gaussSmooth(mesh, tailLeftLight, 1, 0.5);
		
		Vector aLeft = Tools.solveParamInverse(mesh, tailLeftLight, 
				model.getDelta(), model.k,FC.c(0.1));
		Tools.plotVector(mesh, outputFolder, "aLeft"+timeFormat+".dat", aLeft);
		Tools.plotVector(mesh, outputFolder, "aLeftDiff"+timeFormat+".dat", aLeft.add(-1.0, aLeft0));
		//去掉负值，再光滑
		aLeft = Utils.gaussSmooth(mesh, aLeft, 1, 0.5);
		aLeft = Utils.gaussSmooth(mesh, aLeft, 1, 0.4);
		aLeft = Utils.gaussSmooth(mesh, aLeft, 1, 0.3);
		aLeft = Utils.gaussSmooth(mesh, aLeft, 2, 0.4);
		aLeft = Utils.gaussSmooth(mesh, aLeft, 2, 0.3);
		for(int i=1;i<aLeft.getDim();i++) {
			if(aLeft.get(i)<0.0) aLeft.set(i,0.0);
		}
		Tools.plotVector(mesh, outputFolder, "aLeftSmooth_Cut"+timeFormat+".dat", aLeft);

		
		////ampFactor = 1000;
		model.setLightPosition(rightLightCoord[0], rightLightCoord[1]);
		Vector uRight = model.solveNeumann(mesh);
		if(debug) Tools.plotVector(mesh, outputFolder, "uRight"+timeFormat+".dat", uRight);
		
		Vector tailRightLight0 = computeTailRightLight(mesh, uRight);
		if(debug) Tools.plotVector(mesh, outputFolder, "tailRightLight0"+timeFormat+".dat", tailRightLight0);
		Vector aRight0 = Tools.solveParamInverse(mesh, tailRightLight0, 
				model.getDelta(), model.k,FC.c(0.1));
		if(debug) Tools.plotVector(mesh, outputFolder, "aRight0"+timeFormat+".dat", aRight0);
		
		Vector rightRawData = new SparseVectorHashMap(uRight.getDim());
		Vector rightRawDataFull = new SparseVectorHashMap(uRight.getDim());
		for(int i=0;i<rightData.length;i++) {
			Node node = topList.at(i+1);
			if(i>=padRight1 && i<rightData.length-padRight2) {
				uRight.set(node.globalIndex,uRight.get(node.globalIndex)+ampFactor*rightData[i]);
				rightRawData.set(node.globalIndex, rightData[i]);
			}
			rightRawDataFull.set(node.globalIndex, rightData[i]);
		}
		if(debug) Tools.plotVector(mesh, outputFolder, "rightRawData"+timeFormat+".dat", rightRawData);
		if(debug) Tools.plotVector(mesh, outputFolder, "rightRawDataFull"+timeFormat+".dat", rightRawDataFull);

		//Right tail
		Vector tailRightLight = computeTailRightLight(mesh, uRight);
		if(debug) {
			Tools.plotVector(mesh, outputFolder, "tailRightLight"+timeFormat+".dat", tailRightLight);
			Tools.plotVector(mesh, outputFolder, "tailRightLightDiff"+timeFormat+".dat", 
					FMath.axpy(-1, tailRightLight0,tailRightLight));
		}
		//smooth tail
		//tailRightLight = Utils.gaussSmooth(mesh, tailRightLight, 1, 0.5);

		Vector aRight = Tools.solveParamInverse(mesh, tailRightLight, 
				model.getDelta(), model.k,FC.c(0.1));
		Tools.plotVector(mesh, outputFolder, "aRight"+timeFormat+".dat", aRight);		
		Tools.plotVector(mesh, outputFolder, "aRightDiff"+timeFormat+".dat", aRight.add(-1.0, aRight0));
		//去掉负值，再光滑
		aRight = Utils.gaussSmooth(mesh, aRight, 1, 0.5);
		aRight = Utils.gaussSmooth(mesh, aRight, 1, 0.4);
		aRight = Utils.gaussSmooth(mesh, aRight, 1, 0.3);
		aRight = Utils.gaussSmooth(mesh, aRight, 2, 0.4);
		aRight = Utils.gaussSmooth(mesh, aRight, 2, 0.3);
		for(int i=1;i<aRight.getDim();i++) {
			if(aRight.get(i)<0.0) aRight.set(i,0.0);
		}
		Tools.plotVector(mesh, outputFolder, "aRightSmooth_Cut"+timeFormat+".dat", aRight);

		//Scale before average
		double maxLeft = FMath.max(aLeft);
		double minLeft = FMath.min(aLeft);
		double maxRight = FMath.max(aRight);
		double minRight = FMath.min(aRight);
		//按照Right来Scale
		double maxLR = maxRight;//Math.max(maxLeft, maxRight);
		double minLR = minRight;//Math.min(minLeft, minRight);
		//aLeft.scale(maxLR/maxLeft);
		for(int i=1;i<=aLeft.getDim();i++) {
			double t = aLeft.get(i);
			if(t>=0) 
				aLeft.set(i, t*maxLR/maxLeft);
			else
				aLeft.set(i, t*minLR/minLeft);
		}
		//aRight.scale(maxLR/maxRight);
		for(int i=1;i<=aRight.getDim();i++) {
			double t = aRight.get(i);
			if(t>=0) 
				aRight.set(i, t*maxLR/maxRight);
			else
				aRight.set(i, t*minLR/minRight);
		}
		if(debug) {
			Tools.plotVector(mesh, outputFolder, "aLeftSmooth_Cut_Scale"+timeFormat+".dat", aLeft);
			Tools.plotVector(mesh, outputFolder, "aRightSmooth_Cut_Scale"+timeFormat+".dat", aRight);
		}
		Vector aAvg = FMath.axpy(1.0, aLeft, aRight).scale(0.5);
		
		aAvg = Utils.gaussSmooth(mesh, aAvg, 1, 0.5);
		aAvg = Utils.gaussSmooth(mesh, aAvg, 1, 0.4);
		aAvg = Utils.gaussSmooth(mesh, aAvg, 1, 0.3);
		aAvg = Utils.gaussSmooth(mesh, aAvg, 2, 0.4);
		aAvg = Utils.gaussSmooth(mesh, aAvg, 2, 0.3);
		Tools.plotVector(mesh, outputFolder2, String.format("aAvgSmooth_%02d"+timeFormat+".dat",sliceNo), aAvg);
		
//		for(int i=1;i<aAvg.getDim();i++) {
//			if(aAvg.get(i)<0.0) aAvg.set(i,0.0);
//		}
//		Tools.plotVector(mesh, outputFolder2, String.format("aAvgSmoothCut_%02d"+timeFormat+".dat",sliceNo), aAvg);
		
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
	 * T->top data
	 * B->bottom data
	 * @param dataFolder
	 * @return
	 */
	public Map<String,List<double[]>> readAllData(String dataFolder) {
		String[] all = {"T","T"}; //T B
		Map<String,List<double[]>> rlt = new HashMap<String,List<double[]>>();
		for(int j=0;j<all.length;j++) {
			List<double[]> list = new ArrayList<double[]>();
			for(int i=1;i<=NN;i++) {
				String file = String.format(this.outputFolder+"/%s/HumanPhantom2_%s%d.txt", 
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
	 * 3D mesh
	 * [-4.4,4.4]*[-3.8,3.8]*[0,-6]
	 * @param alphaList
	 */
	public void buildResult3D(List<Vector> alphaList,String rltFile) {
		Mesh mesh3D = read3DMesh(this.outputFolder+"/human_phantom3D2.grd");
		Vector v3D = new SparseVectorHashMap(mesh3D.getNodeList().size());
		NodeList nodes = mesh3D.getNodeList();
		double minX = -4.4;
		double dx = 8.8/NN;
		List<Vector2MathFunc> faList = new ArrayList<Vector2MathFunc>();
		for(int i=0;i<alphaList.size();i++) {
			faList.add(new Vector2MathFunc(alphaList.get(i),mesh,"x","y"));
		}
		for(int i=1;i<=nodes.size();i++) {
			Node node = nodes.at(i);
			int index = (int)((node.coord(1) - (minX) - 0.001)/dx);
			v3D.set(i, faList.get(index).apply(
					new Variable("x",node.coord(2)).set("y",node.coord(3))) //3D (x,y,z) -> 2D (y,z)
					);
		}
		Tools.plotVector(mesh3D, outputFolder, rltFile, v3D);
	}
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		HumanPhantom2 hp = new HumanPhantom2();
		hp.init();

		//hp.plot3DInclusion();
		
		//String folder = "phantom-p1";
		String folder = "phantom-p2";
		Map<String,List<double[]>> allData = hp.readAllData(folder);
		
		List<double[]> TD = allData.get("T");
		List<double[]> BD = allData.get("T");
		List<Vector> aList = new ArrayList<Vector>();
		double[] leftLightCoord = {-3.3,-0.5};
		double[] rightLightCoord = {3.3,-0.5};

		//ampFactor需要调整到光强差别在0～10
		for(int i=0;i<TD.size();i++)
			aList.add(hp.run(BD.get(i),TD.get(i),i+1,"out_"+folder,true,
					1000,leftLightCoord,rightLightCoord));
		
		hp.buildResult3D(aList,"3D"+folder+".dat");
	}
	
	public void plot3DInclusion() {
		Mesh mesh3D = read3DMesh("./HumanPhantom/inclusion.grd");
		Vector v3D = new SparseVectorHashMap(mesh3D.getNodeList().size());
		Tools.plotVector(mesh3D, outputFolder, "inclusion.dat", v3D);
	}

}
