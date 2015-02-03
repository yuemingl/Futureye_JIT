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
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.Utils;
import edu.uta.futureye.util.container.NodeList;

public class HumanReal {
	public static String outputFolder = "HumanReal/Output";
	
	String gridFile = "HumanReal/mesh2DSlice.grd";
	int NN = 19; //19*19 grid
	int NSlice = 16;
	Mesh mesh = null;
	NodeList topList = null;
	
	ModelDOT model = new ModelDOT();

	/**
	 * [0,9]*[0,-9]
	 * @param mesh
	 * @param u
	 * @return
	 */
//	public Vector computeTailLeftLight(Mesh mesh, Vector u) {
//		NodeList nodes = mesh.getNodeList();
//		double top = 0.0;
//		double width = 9.0;
//		double skipLen = 1.0; //距离光源1.0范围内的点跳过tail赋值
//		Point lPos = model.getLightPosition();
//		Vector2Function fu = new Vector2Function(u,mesh,"x","y");
//		SparseVector rlt = new SparseVectorHashMap(u.getDim());
//		
//		for(int i=1;i<=nodes.size();i++) {
//			Node node = nodes.at(i);
//			//计算当前点和光源的距离
//			double len = Utils.computeLength(lPos, node);
//			if(len>skipLen && lPos.coord(1)+len < width) {
//				Variable var = new Variable("x",lPos.coord(1)+len).set("y", top);
//				rlt.set(i, fu.value(var));
//			} else {
//				rlt.set(i, u.get(i));
//			}
//		}
//		return rlt;
//	}
	
	/**
	 * [0,9]*[0,-9]
	 * 
	 * 右侧光源在中间
	 * 
	 * @param mesh
	 * @param u
	 * @return
	 */
//	public Vector computeTailRightLight(Mesh mesh, Vector u) {
//		NodeList nodes = mesh.getNodeList();
//		double top = 0;
//		double widthLeft = 6;
//		double widthRight = 3;
//		double skipLen = 1.0;
//		Point rPos = model.getLightPosition();
//		Vector2Function fu = new Vector2Function(u,mesh,"x","y");
//		SparseVector rlt = new SparseVectorHashMap(u.getDim());
//		
//		for(int i=1;i<=nodes.size();i++) {
//			Node node = nodes.at(i);
//			//计算当前点和光源的距离
//			double len = Utils.computeLength(rPos, node);
//			if(node.coord(1) <= rPos.coord(1)) {//光源左侧的点
//				if(len>skipLen && len < widthLeft) {
//					Variable var = new Variable("x",rPos.coord(1)-len).set("y", top);
//					rlt.set(i, fu.value(var));
//				} else {
//					rlt.set(i, u.get(i));
//				}
//			} else {//光源右侧的点
//				if(len>skipLen && len < widthRight) {
//					Variable var = new Variable("x",rPos.coord(1)+len).set("y", top);
//					rlt.set(i, fu.value(var));
//				} else {
//					rlt.set(i, u.get(i));
//				}				
//			}
//			
//		}
//		return rlt;
//	}
	
	public Vector computeTailRightLight(Mesh mesh, Vector u) {
		NodeList nodes = mesh.getNodeList();
		double top = 0.0;
		double xMin = 0.0;
		double xMax = 9.0;
		double skipLen = 1.0;
		Point rPos = model.getLightPosition();
		Vector2Function fu = new Vector2Function(u,mesh,"x","y");
		SparseVector rlt = new SparseVectorHashMap(u.getDim());
		
		for(int i=1;i<=nodes.size();i++) {
			Node node = nodes.at(i);
			//计算当前点和光源的距离
			double len = Utils.computeLength(rPos, node);
			double dx = Math.signum(node.coord(1)-rPos.coord(1))*len*1.0;//!!!
			
			if(len>skipLen && rPos.coord(1)+dx>=xMin && rPos.coord(1)+dx<=xMax) {
				Variable var = new Variable("x",rPos.coord(1)+dx).set("y", top);
				rlt.set(i, fu.apply(var));
			} else {
				rlt.set(i, u.get(i));				
			}
		}
		return rlt;
	}
	public Vector computeTailLeftLight(Mesh mesh, Vector u) {
		return computeTailRightLight(mesh,u);
	}

	public void init() {
		model.setMu_a(FC.c(0.1));
		
		MeshReader reader = new MeshReader(gridFile);
		mesh = reader.read2DMesh();		
		//Use element library to assign degree of freedom (DOF) to element
		Tools.assignLinearShapFunction(mesh);
		mesh.computeNodeBelongsToElements();
		mesh.computeNeighborNodes();
		
		topList = Tools.getNodeListOnY(mesh, 0.0, true);
		//System.out.println(topList);
		
	}

	public Vector run(double[] leftData, double[] rightData, 
			int sliceNo, int timeIndex, String outSubFloder, boolean debug,
			double ampFactor,double[] leftLightCoord, double[] rightLightCoord) {
		
		String outputFolder = String.format("%s/%s/%02d",this.outputFolder,outSubFloder,sliceNo);
		String outputFolder2 = String.format("%s/%s",this.outputFolder,outSubFloder);
		String timeFormat = String.format("_t%05d", timeIndex);

		int padLeft = 6;
		int padRight1 = 7;
		int padRight2 = 4;
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
	
	public Vector run(double[] leftData, double[] rightData, int sliceNo, int timeIndex, String outSubFloder, boolean debug) {
		double[] leftLightCoord = {0.0,0.0};
		double[] rightLightCoord = {6.0,0.0};
		return run(leftData, rightData, 
				sliceNo, timeIndex, outSubFloder, debug,
				1000,leftLightCoord, rightLightCoord);
	}
	
	public Mesh read3DMesh(String file) {
		MeshReader reader = new MeshReader(file);
		Mesh mesh = reader.read3DMesh();
		//Vector v = new SparseVectorHashMap(mesh.getNodeList().size());
		//Tools.plotVector(mesh, outputFolder, "3D.dat", v);
		return mesh;
	}
	
	
	/**
	 * Read all slice of the data
	 * 
	 * Return:
	 * Map of light source and slices list
	 * e.g.
	 * "123" -> slices list (1,...,NSlice)
	 * For each slice we have a sublist of double[], this sublist is a sequence of data in time.
	 * Each slice corresponds to a file, e.g. 
	 * Slice  1: "HumanReal123_1.txt"
	 * Slice  2: "HumanReal123_2.txt"
	 * ...
	 * Slice 16: "HumanReal123_16.txt"
	 * 
	 * For each double[] in the sublists, we have values of light intensity on the points of border
	 * 
	 * @param dataFolder
	 * @return
	 */
	public Map<String,List<List<double[]>>> readAllData(String dataFolder) {
		//Source number
		String[] all = {"123","456","789","101112"};
		
		Map<String,List<List<double[]>>> rlt = new HashMap<String,List<List<double[]>>>();
		for(int j=0;j<all.length;j++) {
			List<List<double[]>> allSlices = new ArrayList<List<double[]>>();
			for(int i=1;i<=NSlice;i++) {
				String file = String.format("./HumanReal/Input/%s/HumanReal%s_%d.txt", 
						dataFolder,	all[j], i);
				List<double[]> oneSlice = readData(file);
				allSlices.add(oneSlice);
			}
			rlt.put(all[j],allSlices);
		}
		return rlt;
	}
	
	public Map<String,List<List<double[]>>> readAllAvgData(String dataFolder) {
		//Source
		String[] all = {"123","456","789","101112"};
		Map<String,List<List<double[]>>> rlt = new HashMap<String,List<List<double[]>>>();
		for(int j=0;j<all.length;j++) {
			List<List<double[]>> allSlices = new ArrayList<List<double[]>>();
			for(int i=1;i<=NSlice;i++) {
				String file = String.format("./HumanReal/Input/%s/HumanRealAvg%s_%d.txt", 
						dataFolder,	all[j], i);
				List<double[]> oneSlice = readData(file);
				allSlices.add(oneSlice);
			}
			rlt.put(all[j],allSlices);
		}
		return rlt;
	}

	/**
	 * Read one slice of the data (one data file)
	 * 
	 * @param fileName
	 * @return list: time sequence;  double[]: values of light intensity on the points of border
	 */
	public List<double[]> readData(String fileName) {
		FileInputStream in;
		try {
			in = new FileInputStream(fileName);

			InputStreamReader reader = new InputStreamReader(in,"UTF-8");
			BufferedReader br = new BufferedReader(reader);
	
			String str = null;
			List<double[]> rlt = new ArrayList<double[]>();
			double[] row=null;
			while((str = br.readLine()) != null){
				String[] line = str.split("(\\s)+");
				row = new double[line.length];
				for(int i=0;i<line.length;i++)
					row[i] = Double.parseDouble(line[i]);
				rlt.add(row);
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
	 * Left:  [0,9]*[0,9]*[ 0, 7.5]
	 * @param alphaSlices
	 * @param rltFile
	 */
	public void buildResult3DLeft(List<Vector> alphaSlices,String rltFile) {
		Mesh mesh3D = read3DMesh("HumanReal/mesh3DLeft.grd");
		Vector v3D = new SparseVectorHashMap(mesh3D.getNodeList().size());
		NodeList nodes = mesh3D.getNodeList();
		double dz = 7.5/(NSlice-1);
		List<Vector2Function> faSlices = new ArrayList<Vector2Function>();
		for(int i=0;i<alphaSlices.size();i++) {
			faSlices.add(new Vector2Function(alphaSlices.get(i),mesh,"x","y"));
		}
		for(int i=1;i<=nodes.size();i++) {
			Node node = nodes.at(i);
			int sliceIndex = (int)((node.coord(3) + 0.001)/dz);
			v3D.set(i, faSlices.get(sliceIndex).apply(
					new Variable("x",node.coord(1)).set("y",node.coord(2)))
					);
		}
		Tools.plotVector(mesh3D, outputFolder, rltFile, v3D);
	}
	
	/**
	 * Right: [0,9]*[0,9]*[10,17.5]
	 * @param alphaSlices
	 * @param rltFile
	 */
	public void buildResult3DRight(List<Vector> alphaSlices,String rltFile) {
		Mesh mesh3D = read3DMesh("HumanReal/mesh3DRight.grd");
		Vector v3D = new SparseVectorHashMap(mesh3D.getNodeList().size());
		NodeList nodes = mesh3D.getNodeList();
		double dz = 7.5/(NSlice-1);
		List<Vector2Function> faSlices = new ArrayList<Vector2Function>();
		for(int i=0;i<alphaSlices.size();i++) {
			faSlices.add(new Vector2Function(alphaSlices.get(i),mesh,"x","y"));
		}
		for(int i=1;i<=nodes.size();i++) {
			Node node = nodes.at(i);
			int sliceIndex = (int)((node.coord(3) -10.0 + 0.001)/dz);
			v3D.set(i, faSlices.get(sliceIndex).apply(
					new Variable("x",node.coord(1)).set("y",node.coord(2)))
					);
		}
		Tools.plotVector(mesh3D, outputFolder, rltFile, v3D);
	}
	
	
	
	public static void printMaxValues(String subFolder,Part part, int T, int prefix) {
		double[] mv = new double[T];
		double[] mi = new double[T];
		String fileName = null;
		for(int t=0;t<T;t++) {
			if(part==Part.LEFT)
				//fileName = String.format(outputFolder+"/"+subFolder+"/3DLeft_t%05d.dat", t);
				fileName = String.format(outputFolder+"/"+subFolder+"/3D_t%05d_Left.dat", t);
			else
				//fileName = String.format(outputFolder+"/"+subFolder+"/3DRight_t%05d.dat", t);
				fileName = String.format(outputFolder+"/"+subFolder+"/3D_t%05d_Right.dat", t);
			Vector v = DataReader.readVector(fileName,4);
			mv[t] = FMath.max(v);
			mi[t] = prefix+t;
		}
		if(part == Part.LEFT)
			DataWriter.writeArray(outputFolder+"/"+subFolder+"/PeakLeft.dat", mi,mv);
		else
			DataWriter.writeArray(outputFolder+"/"+subFolder+"/PeakRight.dat", mi,mv);
			
	}
	
	public static enum Part {LEFT,RIGHT}
	
	/**
	 * Run GCM reconstruction
	 * 
	 * @param dataFolder Output data folder
	 * @param waveFolder Output sub data folder under "dataFolder"
	 * @param part Left or right brain 
	 * @param bAverageData true=do block average
	 */
	public static void run(String dataFolder,String waveFolder, Part part, boolean bAverageData) {
		
		HumanReal hr = new HumanReal();
		hr.init();
		
		Map<String,List<List<double[]>>> allSlices = null;
		//保存Block average的值
		List<List<double[]>> leftSourceBlockAvg = null;
		List<List<double[]>> rightSourceBlockAvg = null;
		
		if(bAverageData) {
			allSlices = hr.readAllAvgData(dataFolder+"/"+waveFolder);
			if(part == Part.LEFT) {
				leftSourceBlockAvg = allSlices.get("123");
				rightSourceBlockAvg = allSlices.get("456");
			} else {
				leftSourceBlockAvg = allSlices.get("789");
				rightSourceBlockAvg = allSlices.get("101112");
			}
		} else {//Do block average here 
			allSlices = hr.readAllData(dataFolder+"/"+waveFolder);
			
			List<List<double[]>> leftSources = null;
			List<List<double[]>> rightSources = null;
			if(part == Part.LEFT) {
				leftSources = allSlices.get("123");
				rightSources = allSlices.get("456");
			} else {
				leftSources = allSlices.get("789");
				rightSources = allSlices.get("101112");
			}
			double dt = 0.092640307993156;
			double[] blockBegin = null;
			if("138".equals(dataFolder)) {
				double[] tmp = {23.4,64.681,101.775,140.057,159.495,212.918,270.48};
				blockBegin = tmp;
			} else if("114".equals(dataFolder)) {
				double[] tmp = {25,67.281,153.453};
				blockBegin = tmp;
			} else if("115".equals(dataFolder)) {
				double[] tmp = {1,48.699,107.808,150.558,173.074,200.089,248.526,293.12,323.713};
				blockBegin = tmp;
			} else if("101".equals(dataFolder)) {
				double[] tmp = {1,44.025,63.025,107.666,152.462,202.352,228.039,333.803};
				blockBegin = tmp;
			} else {
				throw new FutureyeException("Error! Data "+dataFolder+" does not exist!");
			}
	
			double prefix = -5.0;
			double blockLength = 20.0; //200*dt
			//Block average
			leftSourceBlockAvg = new ArrayList<List<double[]>>();
			rightSourceBlockAvg = new ArrayList<List<double[]>>();
			for(int slice=1;slice<=leftSources.size();slice++) {
				List<double[]> left = leftSources.get(slice-1);
				List<double[]> right = rightSources.get(slice-1);
				List<double[]> blockAverageLeft = new ArrayList<double[]>();
				List<Integer> blockAverageLeftCounter = new ArrayList<Integer>();
				List<double[]> blockAverageRight = new ArrayList<double[]>();
				List<Integer> blockAverageRightCounter = new ArrayList<Integer>();
				for(int t=1;t<=left.size();t++) {
					double time = t*dt;
					for(int j=0;j<blockBegin.length;j++) {
						//属于第j个block
						if(time>=blockBegin[j]+prefix && time<blockBegin[j]+blockLength) {
							int index = (int)(time-(blockBegin[j]+prefix));
							
							if(blockAverageLeft.size() < index+1) {
								blockAverageLeft.add(left.get(t));
								blockAverageLeftCounter.add(1);
							} else {
								double[] data = blockAverageLeft.get(index);
								for(int i=0;i<data.length;i++) 
									data[i] += left.get(t)[i];
								blockAverageLeftCounter.set(index, blockAverageLeftCounter.get(index)+1);
							}
							if(blockAverageRight.size() < index+1) {
								blockAverageRight.add(right.get(t));
								blockAverageRightCounter.add(1);
							} else {
								double[] data = blockAverageRight.get(index);
								for(int i=0;i<data.length;i++) 
									data[i] += right.get(t)[i];
								blockAverageRightCounter.set(index, blockAverageRightCounter.get(index)+1);
							}
							break;
						}
					}
				}
				for(int j=0;j<blockAverageLeft.size();j++) {
					double[] data = blockAverageLeft.get(j);
					for(int i=0;i<data.length;i++) 
						data[i] /= blockAverageLeftCounter.get(j);
				}
				for(int j=0;j<blockAverageRight.size();j++) {
					double[] data = blockAverageRight.get(j);
					for(int i=0;i<data.length;i++) 
						data[i] /= blockAverageRightCounter.get(j);
				}
				leftSourceBlockAvg.add(blockAverageLeft);
				rightSourceBlockAvg.add(blockAverageRight);
			}
		}
		
		int time  = leftSourceBlockAvg.get(0).size();
		for(int t=0;t<time;t++) {
			List<Vector> alphaSlices = new ArrayList<Vector>();
			for(int slice=1;slice<=hr.NSlice;slice++) {
				alphaSlices.add(hr.run( leftSourceBlockAvg.get(slice-1).get(t), 
						rightSourceBlockAvg.get(slice-1).get(t),
						slice, t, dataFolder, true));
			}
			if(part==Part.LEFT)
				hr.buildResult3DLeft(alphaSlices, String.format("3D_t%05d_Left.dat",t));
			else
				hr.buildResult3DRight(alphaSlices, String.format("3D_t%05d_Right.dat",t));
		}
		

	}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
//		double[] leftData = new double[19];
//		double[] rightData = new double[19];
//		for(int i=0;i<19;i++) {
//			leftData[i] = -0.01;
//			rightData[i] = -0.01;
//		}
//		
//		HumanReal hr = new HumanReal();
//		hr.init();
//		hr.run(leftData,rightData, 1, 1, "test", true);
		
		//run("138","750avg",Part.LEFT,true);
		//run("138","850avg",Part.LEFT,true);
		//printMaxValues("138-average/750avg",Part.LEFT,25,-4);
		//run("138","750avg",Part.RIGHT,true);
		run("138","850avg",Part.RIGHT,true);
		
		//run("138","750",Part.LEFT);
		//printMaxValues("",Part.LEFT,24,-4);
		//run("138","750",Part.RIGHT);
		//printMaxValues("",Part.RIGHT,24,-4);
		//run("138","850",Part.LEFT);
		//run("138","850",Part.RIGHT);
		
		//run("101","750",Part.LEFT);
		//run("101","750",Part.RIGHT);
		//printMaxValues("101-average",Part.LEFT,24,-4);
		//printMaxValues("101-average",Part.RIGHT,24,-4);
	
		//run("114","750",Part.LEFT);
		//run("114","750",Part.RIGHT);
		//printMaxValues("114-average",Part.LEFT,24,-4);
		//printMaxValues("114-average",Part.RIGHT,24,-4);
		
		//run("115","750",Part.LEFT);
		//run("115","850",Part.LEFT);
		//todo
		//run("115","750",Part.RIGHT);
		//run("115","850",Part.RIGHT);
		//printMaxValues("110-average",Part.LEFT,24,-4);
		//printMaxValues("110-average",Part.RIGHT,24,-4);

	}		

}
