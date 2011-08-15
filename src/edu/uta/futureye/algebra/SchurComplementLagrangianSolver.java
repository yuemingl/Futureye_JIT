package edu.uta.futureye.algebra;

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
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.tutorial.Tools;
import edu.uta.futureye.util.container.DOFList;
import edu.uta.futureye.util.container.ElementList;

/**
 * A = (M  A'  0)
 *     (A  0   C)
 *     (0  C'  R)
 *
 * f = (Fu)
 *     (Fl)
 *     (Fq)
 *
 * x = ( u )
 *     (lmd)
 *     ( p )
 *
 * A' = trans(A)
 * C' = trans(C)
 *
 * 求：A*x=f
 * 
 * (R + C'*inv(A')*M*inv(A)*C)*p 
 *                 = Fq - C'*inv(A')*( Fu - M*inv(A)*Fl )
 *
 * where S = R + C'*inv(A')*M*inv(A)*C, 称为 Schur complement matrix
 *
 *求解Schur补方程，再求解：
 * A*u    = Fl - C*p
 * A'*lmd = Fu - M*u
 * 
 * @author liuyueming
 *
 */
public class SchurComplementLagrangianSolver {
	protected BlockMatrix BM;
	protected BlockVector f;
	protected Mesh mesh;
	
	public SchurComplementLagrangianSolver(BlockMatrix BM,BlockVector f,
			Mesh mesh) {
		this.BM = BM;
		this.f = f;
		this.mesh = mesh;
	}
	
	public BlockVector solve() {
		SparseMatrix M =  (SparseMatrix)BM.getBlock(1, 1);
		SparseMatrix AT = (SparseMatrix)BM.getBlock(1, 2);
		SparseMatrix A =  (SparseMatrix)BM.getBlock(2, 1);
		SparseMatrix C =  (SparseMatrix)BM.getBlock(2, 3);
		SparseMatrix CT = (SparseMatrix)BM.getBlock(3, 2);
		SparseMatrix R =  (SparseMatrix)BM.getBlock(3, 3);
		
		SparseVector Fu = (SparseVector)f.getBlock(1);
		SparseVector Fl = (SparseVector)f.getBlock(2);
		SparseVector Fq = (SparseVector)f.getBlock(3);
		
		FullVector tmp1 = null;
		FullVector tmp2 = null;
		FullVector rhs  = new FullVector(Fu.getDim());
		
		CompressedRowMatrix MM = new CompressedRowMatrix(M, false);
		CompressedRowMatrix RR = new CompressedRowMatrix(R, false);
		CompressedRowMatrix AA = new CompressedRowMatrix(A, false);
		CompressedRowMatrix AAT = new CompressedRowMatrix(AT, false);
		CompressedColMatrix CC = new CompressedColMatrix(C, false);
		CompressedRowMatrix CCT = new CompressedRowMatrix(CT, false);		
		
		FullVector FFu = new FullVector(Fu);
		FullVector FFl = new FullVector(Fl);
		FullVector FFq = new FullVector(Fq);
		
		//Schur complement right hand side: 
		//Fq - C'*inv(A')*( Fu - M*inv(A)*Fl )
		tmp1 = invB_v(AA,FFl);
		MM.mult(tmp1, rhs);
		rhs.axpy(-1.0, FFu);
		
		tmp2 = invB_v(AAT,rhs);
		CCT.mult(tmp2, rhs);
		rhs.axpy(-1.0, FFq);
		
		//S = R + C'*inv(A')*M*inv(A)*C
		//  = R + C'*inv(A*inv(M)*A')*C
		CompressedRowMatrix S = new CompressedRowMatrix();
		CompressedRowMatrix S2 = new CompressedRowMatrix();
		CompressedRowMatrix S3 = new CompressedRowMatrix();
		//R + C'*inv(A')*M*inv(A)*C
		MM.mult(invB_C(AA,CC), S);
		CCT.mult(invB_C(AAT,S.convertToCompressedCol()), S2);
		//R + C'*inv(A*inv(M)*A')*C
		//AA.mult(invB_C(MM,AAT.convertToCompressedCol()),S);
		//S.mult(CC, S3);
		//CCT.mult(S3.convertToCompressedCol(), S2);
		
		S2.axpy(1.0,RR);
		
		Solver sov = new Solver();
		FullVector p = new FullVector(rhs.getDim(),1.0);
		sov.solveCGS(S2, rhs, p);
		//SolverJBLAS sov = new SolverJBLAS();
		//SparseMatrix SS2 = S2.getSparseMatrix();
		//SparseVector Srhs = rhs.getSparseVector();
		//this.imposeDirichletCondition(SS2, Srhs, FC.c0);
		//FullVector p = new FullVector(sov.solveDGESV(SS2, Srhs));
		
		
		//A*u    = Fl - C*p
		//A'*lmd = Fu - M*u
		FullVector u   = new FullVector(Fu.getDim(),1.0);
		FullVector lmd = new FullVector(Fl.getDim(),1.0);
		CC.convertToCompressedRow().mult(p, u);
		u.axpy(-1.0, FFl);
		
		u = invB_v(AA,u);
		//SparseMatrix SAA = AA.getSparseMatrix();
		//SparseVector Su = u.getSparseVector();
		////this.imposeDirichletCondition(SAA, Su, FC.c0);
		//u = new FullVector(sov.solveDGESV(SAA, Su));

		
		MM.mult(u, lmd);
		lmd.axpy(-1.0, FFu);
		lmd = invB_v(AAT,lmd);
		//SparseMatrix SAAT = AAT.getSparseMatrix();
		//SparseVector Sv = v.getSparseVector();
		////this.imposeDirichletCondition(SAAT, Sv, FC.c0);
		//v = new FullVector(sov.solveDGESV(SAAT, Sv));
		
		SparseBlockVector rlt = new SparseBlockVector(3);
		SparseVector uu = new SparseVector(u.getData());
		SparseVector ll = new SparseVector(lmd.getData());
		SparseVector pp = new SparseVector(p.getData());
		rlt.setBlock(1, uu);
		rlt.setBlock(2, ll);
		rlt.setBlock(3, pp);
		return rlt;
	}
	
	
	/**
	 * (M1        |AT1          | 0  )    ( u1 )    ( Fu1 )
	 * (  M2      |   AT2       | 0  )    ( u2 )    ( Fu2 )
	 * (    ...   |      ...    | .. )    ( .. )    ( ..  )
	 * (       MN |         ATN | 0  )    ( uN )    ( FuN )
	 *  ---------------------------        ----       ----
	 * (A1        |0            | C1 )  * (lmd1)  = ( Fl1 )
	 * (  A2      |    0        | C2 )    (lmd2)    ( Fl2 )
	 * (    ...   |      ...    | .. )    ( .. )    ( ..  )
	 * (       AN |          0  | CN )    (lmdN)    ( FlN )
	 *  ---------------------------        ----      ----
	 * (0 0 ... 0 |CT1 CT2 . CTN| R  )    (  q )    ( Fq  )
	 * 
	 * Schur complement:
	 *    S*q = Fq - \sum_{i=1,N}{ CTi*inv(ATi)*( Fui - Mi*inv(Ai)*Fli ) }
	 *  where
	 *    S = R + \sum_{i=1,N}{ CTi*inv(ATi)*Mi*inv(Ai)*Ci }
	 * @return
	 */
	public BlockVector solveMulti() {
		int nBlock = BM.getRowBlockDim();
		int nMeaBlock = (nBlock - 1)/2; //测量次数
		SparseMatrix []M = new SparseMatrix[nMeaBlock];
		SparseMatrix []AT = new SparseMatrix[nMeaBlock];
		SparseMatrix []A = new SparseMatrix[nMeaBlock];
		SparseMatrix []C = new SparseMatrix[nMeaBlock];
		SparseMatrix []CT = new SparseMatrix[nMeaBlock];
		SparseMatrix R = (SparseMatrix)BM.getBlock(nBlock, nBlock);
		SparseVector[] Fu = new SparseVector[nMeaBlock];
		SparseVector[] Fl = new SparseVector[nMeaBlock];
		SparseVector Fq = (SparseVector)f.getBlock(nBlock);		
		
		for(int i=1; i<=nMeaBlock; i++) {
			M[i-1] = (SparseMatrix) BM.getBlock(i, i);
			AT[i-1] = (SparseMatrix) BM.getBlock(i, nMeaBlock+i);
			A[i-1] = (SparseMatrix) BM.getBlock(nMeaBlock+i, i);
			C[i-1] = (SparseMatrix) BM.getBlock(nMeaBlock+i, nBlock);
			CT[i-1] = (SparseMatrix) BM.getBlock(nBlock, nMeaBlock+i);
			
			Fu[i-1] = (SparseVector) f.getBlock(i);
			Fl[i-1] = (SparseVector) f.getBlock(nMeaBlock+i);
		}
		
		CompressedRowMatrix[] MM = new CompressedRowMatrix[nMeaBlock];
		CompressedRowMatrix[] AAT = new CompressedRowMatrix[nMeaBlock];
		CompressedRowMatrix[] AA = new CompressedRowMatrix[nMeaBlock];
		CompressedColMatrix[] CC = new CompressedColMatrix[nMeaBlock];
		CompressedRowMatrix[] CCT = new CompressedRowMatrix[nMeaBlock];
		CompressedRowMatrix RR = new CompressedRowMatrix(R, false);
		
		FullVector[] FFu = new FullVector[nMeaBlock];
		FullVector[] FFl = new FullVector[nMeaBlock];
		FullVector FFq = new FullVector(Fq);
		
		for(int i=0; i<nMeaBlock; i++) {
			MM[i] = new CompressedRowMatrix(M[i], false);
			AAT[i] = new CompressedRowMatrix(AT[i], false);
			AA[i] = new CompressedRowMatrix(A[i], false);
			CC[i] = new CompressedColMatrix(C[i], false);
			CCT[i] = new CompressedRowMatrix(CT[i], false);		
			FFu[i] = new FullVector(Fu[i]);
			FFl[i] = new FullVector(Fl[i]);
		}
		
		//Schur complement right hand side: 
		//Fq - \sum_{i=1,N}{ CTi*inv(ATi)*( Fui - Mi*inv(Ai)*Fli ) }
		FullVector tmp1 = null;
		FullVector tmpSum = new FullVector(Fq.getDim());
		FullVector rhs  = new FullVector(Fq.getDim(),0.0);
		for(int i=0; i<nMeaBlock; i++) {
			tmp1 = invB_v(AA[i],FFl[i]);
			MM[i].mult(tmp1, tmpSum);
			tmpSum.axpy(-1.0, FFu[i]);
		
			tmp1 = invB_v(AAT[i],tmpSum);
			CCT[i].mult(tmp1, tmpSum);
			
			rhs.add(tmpSum);
		}
		rhs.axpy(-1.0, FFq);
		
		
		//S = R + \sum_{i=1,N}{ CTi*inv(ATi)*Mi*inv(Ai)*Ci }
		CompressedRowMatrix tmpS = new CompressedRowMatrix();
		CompressedRowMatrix tmpSSum = new CompressedRowMatrix();
		CompressedRowMatrix S = new CompressedRowMatrix(R.getRowDim(),R.getColDim());
		for(int i=0; i<nMeaBlock; i++) {
			MM[i].mult(invB_C(AA[i],CC[i]), tmpS);
			CCT[i].mult(invB_C(AAT[i],tmpS.convertToCompressedCol()), tmpSSum);
			S.axpy(1.0, tmpSSum);
		}
		S.axpy(1.0,RR);
		//S = RR
		
		Solver sov = new Solver();
		FullVector p = new FullVector(rhs.getDim(),1.0);
		sov.solveCGS(S, rhs, p);
		
		////////////////////////////////////////
		//\Delta(\Delta(sum)) 数值不稳定
		//Vector pp = Tools.computeLaplace2D(mesh, new SparseVector(p.getData()));
		//pp = Tools.computeLaplace2D(mesh, pp);
		//p = new FullVector(pp);
		////////////////////////////////////////
		
		
		//A*ui    = Fli - C*p
		//A'*lmdi = Fui - M*ui
		FullVector[] u = new FullVector[nMeaBlock];
		FullVector[] lmd = new FullVector[nMeaBlock];
		for(int i=0; i<nMeaBlock; i++) {
			u[i] = new FullVector(Fu[i].getDim(),1.0);
			lmd[i] = new FullVector(Fl[i].getDim(),1.0);
			
			CC[i].convertToCompressedRow().mult(p, u[i]);
			u[i].axpy(-1.0, FFl[i]);
			u[i] = invB_v(AA[i],u[i]);
		
			MM[i].mult(u[i], lmd[i]);
			lmd[i].axpy(-1.0, FFu[i]);
			lmd[i] = invB_v(AAT[i],lmd[i]);
		}
		
		SparseBlockVector rlt = new SparseBlockVector(nBlock);
		for(int i=1; i<=nMeaBlock; i++) {
			rlt.setBlock(i, new SparseVector(u[i-1].getData()));
			rlt.setBlock(nMeaBlock+i, new SparseVector(lmd[i-1].getData()));
		}
		rlt.setBlock(nBlock, new SparseVector(p.getData()));
		return rlt;
	}
	
	/**
	 * inv(B)*v
	 * 
	 * @param v
	 * @return
	 */
	public FullVector invB_v(CompressedRowMatrix B, FullVector v) {
		FullVector x = new FullVector(v.getDim(),0.1);
		x.setRandom(1.0, -0.5);
		//SolverJBLAS sov = new SolverJBLAS();
		//x = new FullVector(sov.solveDGESV(B.getSparseMatrix(), v.getSparseVector()));
		Solver sov = new Solver();
		sov.epsRelIter = 1e-5;
		sov.solveCGS(B, v, x);
		return x;
	}
	
	/**
	 * inv(B)*C
	 * @param B
	 * @param C
	 * @return
	 */
	public CompressedColMatrix invB_C(CompressedRowMatrix B, CompressedColMatrix C) {
		CompressedColMatrix BC = new CompressedColMatrix(B.rowDim,C.colDim);
		FullVector v = new FullVector(C.getRowDim());
		FullVector x = new FullVector(C.getRowDim());
		x.setRandom(1.0, -0.5);
		
//		SolverJBLAS sov = new SolverJBLAS();
//		int colDim = C.getColDim();
//		for(int c=1; c<=colDim; c++) {
//			C.getColVector(c, v);
//			x = new FullVector(sov.solveDGESV(B.getSparseMatrix(), v.getSparseVector()));
//			FullVector.SparseData sd = x.getSparseData();
//			BC.setCol(c, sd.index, sd.data);
//		}
		
		Solver sov = new Solver();
		sov.epsRelIter = 1e-5;
		int colDim = C.getColDim();
		for(int c=1; c<=colDim; c++) {
			C.getColVector(c, v);
			sov.solveCGS(B, v, x);
			FullVector.SparseData sd = x.getSparseData();
			BC.setCol(c, sd.index, sd.data);
		}		
		return BC;
	}

	
//	protected void setDirichlet(Matrix BM, Vector BV,
//			int matIndex, double value) {
//		int row = matIndex;
//		int col = matIndex;
//		BM.set(row, col, 1.0);
//		BV.set(row,value);
//		for(int r=1;r<=BM.getRowDim();r++) {
//			if(r != row) {
//				BV.add(r,-BM.get(r, col)*value);
//				BM.set(r, col, 0.0);
//			}
//		}
//		for(int c=1;c<=BM.getColDim();c++) {
//			if(c != col) {
//				BM.set(row, c, 0.0);
//			}
//		}
//	}
//	
//	public void imposeDirichletCondition(SparseMatrix BM, SparseVector BV,
//			Function diri) {
//		ElementList eList = mesh.getElementList();
//		int nNode = mesh.getNodeList().size();
//		for(int i=1;i<=eList.size();i++) {
//			Element e = eList.at(i);
//			DOFList DOFs = e.getAllDOFList(DOFOrder.NEFV);
//			for(int j=1;j<=DOFs.size();j++) {
//				DOF dof = DOFs.at(j);
//				GeoEntity ge = dof.getOwner();
//				if(ge instanceof Node) {
//					Node n = (Node)ge;
//					if(n.getNodeType() == NodeType.Dirichlet) {
//						Variable v = Variable.createFrom(diri, n, 0);
//						setDirichlet(BM,BV,dof.getGlobalIndex(),diri.value(v));
//					}
//				}
//			}
//		}
//	}
}
