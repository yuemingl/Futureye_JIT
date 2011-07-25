package edu.uta.futureye.application;

import java.io.File;
import java.util.HashMap;

import edu.uta.futureye.algebra.SchurComplementLagrangianSolver;
import edu.uta.futureye.algebra.Solver;
import edu.uta.futureye.algebra.SolverJBLAS;
import edu.uta.futureye.algebra.SparseBlockMatrix;
import edu.uta.futureye.algebra.SparseBlockVector;
import edu.uta.futureye.algebra.SparseMatrix;
import edu.uta.futureye.algebra.SparseVector;
import edu.uta.futureye.algebra.intf.BlockMatrix;
import edu.uta.futureye.algebra.intf.BlockVector;
import edu.uta.futureye.algebra.intf.Matrix;
import edu.uta.futureye.algebra.intf.Vector;
import edu.uta.futureye.core.DOF;
import edu.uta.futureye.core.DOFOrder;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.core.NodeType;
import edu.uta.futureye.core.geometry.GeoEntity;
import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.basic.Vector2Function;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.operator.FMath;
import edu.uta.futureye.io.MeshReader;
import edu.uta.futureye.io.MeshWriter;
import edu.uta.futureye.lib.assembler.AssemblerScalar;
import edu.uta.futureye.lib.element.FELinearTriangle;
import edu.uta.futureye.lib.weakform.WeakFormL22D;
import edu.uta.futureye.lib.weakform.WeakFormLaplace2D;
import edu.uta.futureye.util.container.DOFList;
import edu.uta.futureye.util.container.ElementList;
import edu.uta.futureye.util.container.NodeList;

/**
 * Implementation of the paper: 
 *   'A Framework For The Adaptive Finite Element Solution Of Large-Scale Inverse Problems'
 * 
 * Lagrange Multiplier Method
 * 
 * @author liuyueming
 *
 */
public class VariationGaussNewton {
	public Mesh mesh;
	protected static String outputFolder = "Lagrangian_Wolfgang";
	public boolean debug = false;

    Function coef_q = new AbstractFunction("x","y") {
    	@Override
    	public double value(Variable v) {
    		double x = v.get("x");
    		double y = v.get("y");
    		if(Math.sqrt(x*x+y*y)<0.5)
    			return 1.0;
    		else
    			return 8.0;
    	}
    };

    Function diri = new AbstractFunction("x","y") {
    	@Override
    	public double value(Variable v) {
    		double x = v.get("x");
    		double y = v.get("y");
    		return (x*x+y*y)/8.0 + 7.0/32.0;
    	}
    };
    
    //测量数据
    Function z = null;
    Function qBar = null;
    
    //正则化参数
    //double beta = 0.01;
  //  double beta = 0.5;
    double beta = 0.001;
    
    private int iterNum = 0;
    
	public static void plotVector(Mesh mesh, Vector v, String fileName) {
	    MeshWriter writer = new MeshWriter(mesh);
	    if(!outputFolder.isEmpty()) {
		    File file = new File("./"+outputFolder);
			if(!file.exists()) {
				file.mkdirs();
			}
	    }
	    writer.writeTechplot("./"+outputFolder+"/"+fileName, v);
	}

	public static void plotFunction(Mesh mesh, Function fun, String fileName) {
	    NodeList list = mesh.getNodeList();
	    int nNode = list.size();
		Variable var = new Variable();
		Vector v = new SparseVector(nNode);
	    for(int i=1;i<=nNode;i++) {
	    	Node node = list.at(i);
	    	var.setIndex(node.globalIndex);
	    	var.set("x", node.coord(1));
	    	var.set("y", node.coord(2));
	    	v.set(i, fun.value(var));
	    }
	    plotVector(mesh,v,fileName);
	}	
	

	public void readMesh(){
        MeshReader reader = new MeshReader("lagrangian.grd");
        mesh = reader.read2DMesh();
        mesh.computeNodeBelongsToElements();

        //2.Mark border types
        HashMap<NodeType, Function> mapNTF =
                new HashMap<NodeType, Function>();
        mapNTF.put(NodeType.Dirichlet, null);
        mesh.markBorderNode(mapNTF);
        ElementList eList = mesh.getElementList();
//        for(int i=1;i<=eList.size();i++) {
//        	System.out.println(eList.at(i));
//        }

        //3.Use element library to assign degrees of
        //  freedom (DOF) to element
        FELinearTriangle feLT = new FELinearTriangle();
        for(int i=1;i<=eList.size();i++)
            feLT.assignTo(eList.at(i));
        
	}
	
	/**
	 * State equation
	 * L_{\lambda} := (q\nabla{u},\nabla{\phi})-(f,\phi)=0
	 * 
	 * 获取状态方程的余量
	 * @return
	 */
	public Vector getResLlmd(Vector u, Vector q) {
        //4.Weak form
        WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
        Function fq = new Vector2Function(q);
        weakForm.setParam(fq, FC.c0, null, null);
        //Right hand side(RHS): f(x) = -4.0
        weakForm.setF(FC.c(-4.0));

        //5.Assembly process
        AssemblerScalar assembler =
                new AssemblerScalar(mesh, weakForm);
        System.out.println("Begin Assemble...");
        assembler.assemble();
        Matrix stiff = assembler.getStiffnessMatrix();
        Vector load = assembler.getLoadVector();
        //Boundary condition
        assembler.imposeDirichletCondition(diri);
        System.out.println("Assemble done!");

        if(debug) {
	        //6.Solve linear system
	        Solver solver = new Solver();
	        Vector u_solve = solver.solveCGS(stiff, load);
	        plotVector(mesh,u_solve,"u_solve.dat");
        }
        
        //Residual
        Vector res = new SparseVector(u.getDim());
        stiff.mult(u, res);
        res.add(-1.0, load);
        plotVector(mesh,res,"Res_Llambda.dat");
        
        return res;
	}
	
	/**
	 * Adjoint equation
	 * L_{u} := (q\nabla{\psi},\nabla{\lambda}) + (u-z,\psi)
	 * 
	 * @return
	 */
	public Vector getResLu(Vector u, Vector lambda, Vector q) {
        //4.Weak form
        WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
        Function fq = new Vector2Function(q);
        weakForm.setParam(fq, FC.c0, null, null);
        //Right hand side(RHS): f(x) = - (u - z)
        Function z_u = z.S(new Vector2Function(u));
        plotFunction(mesh, z_u, String.format("z_u%02d.dat",this.iterNum));
        weakForm.setF(z_u);

        //5.Assembly process
        AssemblerScalar assembler =
                new AssemblerScalar(mesh, weakForm);
        System.out.println("Begin Assemble...");
        assembler.assemble();
        Matrix stiff = assembler.getStiffnessMatrix();
        Vector load = assembler.getLoadVector();
        //Boundary condition
        assembler.imposeDirichletCondition(FC.c0);
        System.out.println("Assemble done!");

        if(debug) {
	        //6.Solve linear system
	        Solver solver = new Solver();
	        Vector lmd_solve = solver.solveCGS(stiff, load);
	        plotVector(mesh, lmd_solve, "lambda_solve.dat");
        }
        
        //Residual
        Vector res = new SparseVector(lambda.getDim());
        stiff.mult(lambda, res);
        res.add(-1.0, load);
        plotVector(mesh,res,"Res_Lu.dat");
        
        return res;
    }
	
	/**
	 * Parameter regularization
	 * L_{q} := \beta(q-\overline{q},\chi) + (\chi\nabla{u},\nabla{\lambda})
	 * 
	 * @param u
	 * @param lambda
	 * @param q
	 * @return
	 */
	public Vector getResLq(Vector u, Vector lambda, Vector q) {
        //4.Weak form
        WeakFormL22D weakForm = new WeakFormL22D();
        weakForm.setParam(FC.c0, FC.c1);
        //Right hand side(RHS): f(x) = -(1.0/\beta)\nabla{u}\cdot\nabla{v}
        Function fu = new Vector2Function(u,mesh,"x","y");
        Function flmd = new Vector2Function(lambda,mesh,"x","y");
        Function f = FMath.grad(fu).dot(FMath.grad(flmd));
        plotFunction(mesh,f,String.format("Grad(u)Grad(lmd)%02d.dat",this.iterNum));
        Function f2 = FC.c(-1.0/beta).M(f).A(qBar);
        plotFunction(mesh,f2,String.format("LqRHS%02d.dat",this.iterNum));
        weakForm.setF(f2);

        //5.Assembly process
        AssemblerScalar assembler =
                new AssemblerScalar(mesh, weakForm);
        System.out.println("Begin Assemble...");
        assembler.assemble();
        Matrix stiff = assembler.getStiffnessMatrix();
        Vector load = assembler.getLoadVector();
        //Boundary condition
        assembler.imposeDirichletCondition(FC.c(8.0));
        System.out.println("Assemble done!");

        if(debug) {
	        //6.Solve linear system
	        Solver solver = new Solver();
	        Vector q_solve = solver.solveCGS(stiff, load);
	        plotVector(mesh, q_solve, "q_solve.dat");
        }
        
        //Residual
        Vector res = new SparseVector(q.getDim());
        stiff.mult(q, res);
        res.add(-1.0, load);
        plotVector(mesh,res,"Res_Lq.dat");
        
        return res;
	}
	
	
	
	/**
	 * Search direction matrix:
	 * 
	 *  ( M  A'  0 )(du)     (Lu)
	 *  ( A  0   C )(dl) = - (Ll)
	 *  ( 0  C' bR )(dq)     (Lq)
	 * 
	 * @return
	 */
	public Matrix getM() {
        //4.Weak form: (du,\phi)
        WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
        weakForm.setParam(FC.c0, FC.c1, null, null);
        //Right hand side(RHS): f(x) = 0
        weakForm.setF(FC.c0);

        //5.Assembly process
        AssemblerScalar assembler =
                new AssemblerScalar(mesh, weakForm);
        System.out.println("Begin Assemble...");
        assembler.assemble();
        Matrix stiff = assembler.getStiffnessMatrix();
        //Boundary condition 
        //边界条件需要在大矩阵中设置
        //assembler.imposeDirichletCondition(FC.c0);
        System.out.println("Assemble done!");
        
        return stiff;
    }
	
	public Matrix getA(Vector qk) {
        //4.Weak form: (\nabla{dl},qk*\nabla{\phi})
        WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
        Function fqk = new Vector2Function(qk);
        weakForm.setParam(fqk, FC.c0, null, null);
        //Right hand side(RHS): f(x) = 0
        weakForm.setF(FC.c0);

        //5.Assembly process
        AssemblerScalar assembler =
                new AssemblerScalar(mesh, weakForm);
        System.out.println("Begin Assemble...");
        assembler.assemble();
        Matrix stiff = assembler.getStiffnessMatrix();
        //Boundary condition
        //边界条件需要在大矩阵中设置
        //assembler.imposeDirichletCondition(FC.c0);
        System.out.println("Assemble done!");
        
        return stiff;
	}
	
	public Matrix getC(Vector uk) {
        //4.Weak form: (\nabla{\phi},dq\nabla{uk})
        WeakFormGCM weakForm = new WeakFormGCM();
        
        Function fuk = new Vector2Function(uk,mesh,"x","y");
        Function u_x = fuk._d("x");
        Function u_y = fuk._d("y");
        plotFunction(mesh, u_x, "u_x.dat");
        plotFunction(mesh, u_y, "u_y.dat");
        weakForm.setParam(FC.c0, FC.c0, u_x, u_y);
        //Right hand side(RHS): f(x) = 0
        weakForm.setF(FC.c0);

        //5.Assembly process
        AssemblerScalar assembler =
                new AssemblerScalar(mesh, weakForm);
        System.out.println("Begin Assemble...");
        assembler.assemble();
        Matrix stiff = assembler.getStiffnessMatrix();
        //Boundary condition
        //边界条件需要在大矩阵中设置
        //assembler.imposeDirichletCondition(FC.c0);
        System.out.println("Assemble done!");
        
        return stiff;
	}
	
	public Matrix getBR() {
        //4.Weak form: (dq,\chi)
        WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
        weakForm.setParam(FC.c0, FC.c(beta), null, null);
        //Right hand side(RHS): f(x) = 0
        weakForm.setF(FC.c0);

        //5.Assembly process
        AssemblerScalar assembler =
                new AssemblerScalar(mesh, weakForm);
        System.out.println("Begin Assemble...");
        assembler.assemble();
        Matrix stiff = assembler.getStiffnessMatrix();
        //Boundary condition
        //边界条件需要在大矩阵中设置
        //assembler.imposeDirichletCondition(FC.c0);
        System.out.println("Assemble done!");
        
        return stiff;
	}
	
	public Vector getRealU() {
        //4.Weak form
        WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
        plotFunction(mesh,coef_q,"qReal.dat");
        weakForm.setParam(coef_q, FC.c0, null, null);
        //Right hand side(RHS): f(x) = -4.0
        weakForm.setF(FC.c(-4.0));

        //5.Assembly process
        AssemblerScalar assembler =
                new AssemblerScalar(mesh, weakForm);
        System.out.println("Begin Assemble...");
        assembler.assemble();
        Matrix stiff = assembler.getStiffnessMatrix();
        Vector load = assembler.getLoadVector();
        //Boundary condition
        assembler.imposeDirichletCondition(diri);
        System.out.println("Assemble done!");

        //6.Solve linear system
        Solver solver = new Solver();
        Vector u = solver.solveCG(stiff, load);
        System.out.println("u=");
        for(int i=1;i<=u.getDim();i++)
            System.out.println(String.format("%.3f", u.get(i)));

        plotVector(mesh,u,"uReal.dat");
        return u;
	}
	
	public Vector solveStateEquation(Vector q) {
        //4.Weak form
        WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
        Function fq = new Vector2Function(q);
        weakForm.setParam(fq, FC.c0, null, null);
        //Right hand side(RHS): f(x) = -4.0
        weakForm.setF(FC.c(-4.0));

        //5.Assembly process
        AssemblerScalar assembler =
                new AssemblerScalar(mesh, weakForm);
        System.out.println("Begin Assemble...");
        assembler.assemble();
        Matrix stiff = assembler.getStiffnessMatrix();
        Vector load = assembler.getLoadVector();
        
        //Boundary condition
        assembler.imposeDirichletCondition(diri);//???给定一个初始猜测q,边界条件是什么？
        System.out.println("Assemble done!");

        //6.Solve linear system
        Solver solver = new Solver();
        Vector u_solve = solver.solveCGS(stiff, load);
       
        return u_solve;
	}
	
	public Vector solveAdjointEquation(Vector u, Vector q) {
        //4.Weak form
        WeakFormLaplace2D weakForm = new WeakFormLaplace2D();
        Function fq = new Vector2Function(q);
        weakForm.setParam(fq, FC.c0, null, null);
        //Right hand side(RHS): f(x) = - (u - z)
        weakForm.setF(z.S(new Vector2Function(u)));

        //5.Assembly process
        AssemblerScalar assembler =
                new AssemblerScalar(mesh, weakForm);
        System.out.println("Begin Assemble...");
        assembler.assemble();
        Matrix stiff = assembler.getStiffnessMatrix();
        Vector load = assembler.getLoadVector();
        //Boundary condition
        assembler.imposeDirichletCondition(FC.c0);
        System.out.println("Assemble done!");

        //6.Solve linear system
        Solver solver = new Solver();
        Vector lmd_solve = solver.solveCGS(stiff, load);
        
        return lmd_solve;
    }
	
	protected void setDirichlet(BlockMatrix BM, BlockVector BV,
			int matIndex, double value) {
		int row = matIndex;
		int col = matIndex;
		BM.set(row, col, 1.0);
		BV.set(row,value);
		for(int r=1;r<=BM.getRowDim();r++) {
			if(r != row) {
				BV.add(r,-BM.get(r, col)*value);
				BM.set(r, col, 0.0);
			}
		}
		for(int c=1;c<=BM.getColDim();c++) {
			if(c != col) {
				BM.set(row, c, 0.0);
			}
		}
	}
	
	public void imposeDirichletCondition(BlockMatrix BM, BlockVector BV,
			Function diri) {
		ElementList eList = mesh.getElementList();
		int nNode = mesh.getNodeList().size();
		for(int i=1;i<=eList.size();i++) {
			Element e = eList.at(i);
			DOFList DOFs = e.getAllDOFList(DOFOrder.NEFV);
			for(int j=1;j<=DOFs.size();j++) {
				DOF dof = DOFs.at(j);
				GeoEntity ge = dof.getOwner();
				if(ge instanceof Node) {
					Node n = (Node)ge;
					if(n.getNodeType() == NodeType.Dirichlet) {
						Variable v = Variable.createFrom(diri, n, 0);
						setDirichlet(BM,BV,dof.getGlobalIndex(),diri.value(v));
						setDirichlet(BM,BV,nNode+dof.getGlobalIndex(),diri.value(v));
						setDirichlet(BM,BV,nNode*2+dof.getGlobalIndex(),diri.value(v));
					}
				}
			}
		}
	}
	
	public BlockMatrix getSearchDirectionMatrix(Vector uk, Vector qk) {
		BlockMatrix BM = new SparseBlockMatrix(3,3);
		Matrix M = this.getM();
		Matrix A = this.getA(qk);
		Matrix AT = A.copy().trans();
		Matrix C = this.getC(uk);
		Matrix CT = C.copy().trans();
		Matrix R = this.getBR();
		
		Matrix BM13 = new SparseMatrix(M.getRowDim(),R.getColDim());
		Matrix BM22 = new SparseMatrix(A.getRowDim(),A.getColDim());
		Matrix BM31 = new SparseMatrix(R.getRowDim(),M.getColDim());
		
		BM.setBlock(1, 1, M);
		BM.setBlock(1, 2, AT);
		BM.setBlock(1, 3, BM13);
		
		BM.setBlock(2, 1, A);
		BM.setBlock(2, 2, BM22);
		BM.setBlock(2, 3, C);
		
		BM.setBlock(3, 1, BM31);
		BM.setBlock(3, 2, CT);
		BM.setBlock(3, 3, R);
		
		return BM;
	}
	
	/**
	 * Return a new BolckMatrix object that share the same sub matrix objects
	 * 
	 * @param BM
	 * @param col1
	 * @param col2
	 * @return
	 */
	public BlockMatrix changeBlockColumn(BlockMatrix BM, int col1, int col2) {
		int blockRow = BM.getRowBlockDim();
		int blockCol = BM.getColBlockDim();
		
		BlockMatrix newBM = new SparseBlockMatrix(blockRow, blockCol);
		for(int i=1;i<=blockRow;i++) {
			for(int j=1;j<=blockCol;j++) {
				newBM.setBlock(i, j, BM.getBlock(i, j));
			}
		}
		
		for(int i=1;i<=BM.getRowBlockDim();i++) {
			newBM.setBlock(i, col1, BM.getBlock(i, col2));
			newBM.setBlock(i, col2, BM.getBlock(i, col1));
		}
		return newBM;
	}	
	
	
	public void gaussNewtonIterate() {
		//初值
		NodeList nodes = mesh.getNodeList();
		//初始q=q0
		Vector q0 = new SparseVector(nodes.size());
		Vector q00 = new SparseVector(nodes.size());
		for(int i=1;i<=nodes.size();i++) {
			Node node = nodes.at(i);
			if(node.getNodeType()==NodeType.Dirichlet)
				q0.set(i,this.coef_q.value(Variable.createFrom(coef_q, node, 0)));
			else {
				//q0.set(i,this.coef_q.value(Variable.createFrom(coef_q, node, 0))+1.0+Math.random()*0.9);
				//q0.set(i,1.0+Math.random()*0.9);
				
				double x = node.coord(1);
				double y = node.coord(2);
				
				if(Math.sqrt((x)*(x)+(y)*(y))<0.5) {
					q0.set(i,3.0);
				
				//=1的圆位置偏移
				//if(Math.sqrt((x-0.3)*(x-0.3)+(y)*(y))<0.5) {
				//	q0.set(i,1.0);
				} else {
					q0.set(i,8.0);
				}
			}
			q00.set(i,this.coef_q.value(Variable.createFrom(coef_q, node, 0))+1.0+Math.random()*0.9);
				
		}
		plotVector(mesh, q0, "q0.dat");
		plotVector(mesh, q00, "q00.dat");
		
		//初始u=u0, lambda=lambda0，从q00计算而来
		Vector u0 = this.solveStateEquation(q0);//q00
		plotVector(mesh, u0, "u0.dat");
		Vector lambda0 = this.solveAdjointEquation(u0, q0);//q00
		plotVector(mesh, lambda0, "lambda0.dat");
		
		//强制u0=0,lambda0=0
		//q0 = new SparseVector(nodes.size(),1.0);
		//u0 = new SparseVector(nodes.size(),0.0);
		//lambda0 = new SparseVector(nodes.size(),0.0);
		
		
		BlockVector f = new SparseBlockVector(3);
		//迭代求解
		for(int iter=0;iter<50;iter++) {
			iterNum = iter;
			
			BlockMatrix BM = this.getSearchDirectionMatrix(u0,q0);
			
			f.setBlock(1, this.getResLu(u0, lambda0, q0).scale(-1.0));
			f.setBlock(2, this.getResLlmd(u0, q0).scale(-1.0));//bugfix lambda0->q0
			f.setBlock(3, this.getResLq(u0, lambda0, q0).scale(-1.0));

			//设置边界条件并求解
			//需要交换矩阵BM的第1, 2列，然后设置边界条件，这样使得矩阵A, A'可逆
			//BlockMatrix newBM = this.changeBlockColumn(BM, 1, 2);
			//this.imposeDirichletCondition(newBM, f, FC.c0);
			//SchurComplementLagrangianSolver solver = 
			//	new SchurComplementLagrangianSolver(BM, f,mesh);
			//BlockVector x = solver.solve();
			
			this.imposeDirichletCondition(BM, f, FC.c0);
			SolverJBLAS sol = new SolverJBLAS();
			BlockVector x = (BlockVector)sol.solveDGESV(BM, f);
			
			plotVector(mesh, x.getBlock(1), String.format("delta_u%02d.dat",iter));
			plotVector(mesh, x.getBlock(2), String.format("delta_lambda%02d.dat",iter));
			plotVector(mesh, x.getBlock(3), String.format("delta_q%02d.dat",iter));
			
			//待定，与beta选取有关
			//double stepLength = 0.01;
			//double stepLength = 0.3;
			//double stepLength = 0.00005;
			double stepLength = 0.0001;
			
			
			u0.add(stepLength, x.getBlock(1));
			lambda0.add(stepLength, x.getBlock(2));
			q0.add(stepLength, x.getBlock(3));
			
			plotVector(mesh, u0, String.format("rlt_u%02d.dat",iter));
			plotVector(mesh, lambda0, String.format("rlt_lambda%02d.dat",iter));
			plotVector(mesh, q0, String.format("rlt_q%02d.dat",iter));

		}

		
	}
	
	public static void main(String[] args) {
		VariationGaussNewton vgn = new VariationGaussNewton();
		vgn.debug = true;
		vgn.readMesh();
		
		
		Vector uReal = vgn.getRealU();
//		NodeList nodes = vgn.mesh.getNodeList();
//		for(int i=1;i<=nodes.size();i++) {
//			Node node = nodes.at(i);
//			if(node.coord(2)<0.0)
//				uReal.set(i,0.0);
//		}
//		plotVector(vgn.mesh, uReal, "uRealCut.dat");
		
		//u测量值
		vgn.z = new Vector2Function(uReal);
		//q参考值
		vgn.qBar = new AbstractFunction("x","y") {
	    	@Override
	    	public double value(Variable v) {
	    		double x = v.get("x");
	    		double y = v.get("y");
	    		
	    		//如果参考值=1.0的圆的值增加到3.0
	    		if(Math.sqrt(x*x+y*y)<0.5)
	    			return 3.0;
	    		
	    		
	    		//如果参考值=1.0的圆位置偏移，结果怎么样？
	    		//结果会偏向参考值
		    	//if(Math.sqrt((x-0.3)*(x-0.3)+y*y)<0.5)
    			//	return 1.0;
	    		
	    		//如果参考值=1.0的圆的中间更小的圆的值增加到3.0
	    		//if(Math.sqrt(x*x+y*y)<0.3) {
	    		//	return 3.0;
	    		//} else if(Math.sqrt(x*x+y*y)<0.5)
    			//	return 1.0;
	    		
	    		else
	    			return 8.0;
	    	}
	    };
	    VariationGaussNewton.plotFunction(vgn.mesh, vgn.qBar, "qBar.dat");
		
/*Test	    
		Vector u0 = new SparseVector(vgn.mesh.getNodeList().size(),1.0);
		//u0 = uReal;
		Vector lambda0 = new SparseVector(vgn.mesh.getNodeList().size(),1.0);
		Vector q0 = new SparseVector(vgn.mesh.getNodeList().size(),1.0);
		Vector resLlmd = vgn.getResLlmd(u0, q0);
		Vector resLu = vgn.getResLu(u0, lambda0, q0);
		Vector resLq = vgn.getResLq(u0, lambda0, q0);
*/		
		//bugfix 修改q的边界条件1.0->8.0
		//bugfix (Lu,Ll,Lq)右端项有负号
		//bugfix getC() qk->0
		//bugfix getRes*() res.axpy(-1.0, load); -> res.add(-1.0, load);
	    //bugfix Solver.solveCGS()默认不要清空系数矩阵
		vgn.gaussNewtonIterate();
		
	}
}
