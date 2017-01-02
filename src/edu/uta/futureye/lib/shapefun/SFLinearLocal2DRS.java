package edu.uta.futureye.lib.shapefun;

import static edu.uta.futureye.function.FMath.C1;
import static edu.uta.futureye.function.FMath.r;
import static edu.uta.futureye.function.FMath.s;

import java.util.HashMap;
import java.util.Map;

import org.apache.bcel.generic.ConstantPoolGen;
import org.apache.bcel.generic.InstructionFactory;
import org.apache.bcel.generic.InstructionHandle;
import org.apache.bcel.generic.InstructionList;
import org.apache.bcel.generic.MethodGen;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.function.MultiVarFunc;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.function.intf.ScalarShapeFunction;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.Utils;
import edu.uta.futureye.util.container.ObjList;
import edu.uta.futureye.util.container.VertexList;

/**
 * Linear shape function in triangle local coordinate system (area coordinate system)
 * 
 *   Ni = Ni(r,s) = Ni( r(x,y), s(x,y) ), i=1,2,3
 *   
 * where
 *     N1 = r
 *     N2 = s
 *     N3 = 1 - r - s
 */
public class SFLinearLocal2DRS  extends MultiVarFunc 
							  implements ScalarShapeFunction {
	protected int funIndex;
	private MathFunc funCompose = null;
	private MathFunc funOuter = null;
	protected ObjList<String> innerVarNames = null;
	
	protected Element e = null;
	private double jac = 0.0;
	private double[] x = new double[3];
	private double[] y = new double[3];
	
	private double coef = 1.0;
	
	/**
	 * 构造下列形函数中的一个：
	 * N1 = L1 = r
	 * N2 = L2 = s
	 * N3 = L3 = t = 1 - r - s
	 * @param funID = 1,2,3
	 * 
	 */
	public void Create(int funID, double coef) {
		funIndex = funID - 1;
		if( funID<1 || 3<funID ) {
			throw new FutureyeException("ERROR: funID should be 1,2 or 3.");
		}
		
		this.varNames = new String[]{"r","s"};
		innerVarNames = new ObjList<String>("x","y");
		
		//Compose function: r=r(x,y), s=s(x,y)
		Map<String, MathFunc> fInners = new HashMap<String, MathFunc>();
		for(final String varName : varNames) {
			fInners.put(varName, new MultiVarFunc(varName, innerVarNames.toList()) {	
				public MathFunc diff(String var) {
					if(varName.equals("r")) { //r对应三角形高h的负倒数
						if(var.equals("x"))
							return FC.c( (y[1]-y[2]) / jac);
						if(var.equals("y"))
							return FC.c( (x[2]-x[1]) / jac);
					} else if(varName.equals("s")) { //s对应三角形高h的负倒数
						if(var.equals("x"))
							return FC.c( (y[2]-y[0]) / jac);
						if(var.equals("y"))
							return FC.c( (x[0]-x[2]) / jac);
					} else
						throw new FutureyeException("\nERROR:\n varName="+varName);
					return null;
				}
				@Override
				public double apply(double... args) {
					throw new UnsupportedOperationException();
				}
				
				@Override
				public String getExpr() {
					return varName;
				}
				@Override
				public String toString() {
					return varName+"(x,y)";
				}
				
			});
		}
		
		//Construct shape functions: r,s are free variables, t = 1 - r - s
		if(funIndex == 0) {
			funOuter = r; //N1 = r
		} else if(funIndex == 1) {
			funOuter = s; //N2 = s
		} else { 
			funOuter = C1.S(r).S(s); //N3 = t = 1 - r - s
		}
		this.coef = coef;
		funCompose = funOuter.compose(fInners).M(FC.c(this.coef));
		/**
		 * The default active variable names of a composite function is the inner variable names.
		 * Shape function needs the outer variable names as the active variable names.
		 */
		funCompose.setActiveVarNames(funOuter.getVarNames());
		funCompose.setArgIdx(Utils.getIndexMap(this.getVarNames()));
	}
	
	public SFLinearLocal2DRS(int funID) {
		this.Create(funID, 1.0);
	}
	
	public SFLinearLocal2DRS(int funID, double coef) {
		this.Create(funID, coef);
	}
	
	@Override
	public MathFunc diff(String varName) {
		if(e == null) {
			throw new FutureyeException("Call assignElement first before calling diff(\""+varName+"\")!");
		}
		return funCompose.diff(varName);
	}

	@Override
	public double apply(Variable v) {
		return funCompose.apply(v);
	}
	
	@Override
	public double apply(Element e, Node n, double ...args) {
		return funCompose.apply(e, n, args);
	}
	
	@Override
	public InstructionHandle bytecodeGen(String clsName, MethodGen mg,
			ConstantPoolGen cp, InstructionFactory factory,
			InstructionList il, Map<String, Integer> argsMap, int argsStartPos, 
			Map<MathFunc, Integer> funcRefsMap) {
		return funCompose.bytecodeGen(clsName, mg, cp, factory, il, argsMap, argsStartPos, funcRefsMap);
	}
	
	@Override
	public MathFunc setArgIdx(Map<String, Integer> argsMap) {
		super.setArgIdx(argsMap);
		this.funCompose.setArgIdx(argsMap);
		return this;
	}
	
	@Override
	public void assignElement(Element e) {
		this.e = e;
		VertexList vList = e.vertices();
		x[0] = vList.at(1).coord(1);
		x[1] = vList.at(2).coord(1);
		x[2] = vList.at(3).coord(1);
		y[0] = vList.at(1).coord(2);
		y[1] = vList.at(2).coord(2);
		y[2] = vList.at(3).coord(2);
		jac = (x[0]-x[2])*(y[1]-y[2])-(x[1]-x[2])*(y[0]-y[2]);
	}

	public String getExpr() {
		return "N"+(funIndex+1)+"(r,s)";
	}
	
	public String toString() {
		return "N"+(funIndex+1)+"(r,s) = "+funOuter.getExpr();
	}
	
	ScalarShapeFunction sf1d1 = new SFLinearLocal1D(1);
	ScalarShapeFunction sf1d2 = new SFLinearLocal1D(2);
	@Override
	public ScalarShapeFunction restrictTo(int funIndex) {
		if(funIndex == 1) return sf1d1;
		else return sf1d2;
	}

	@Override
	public ObjList<String> innerVarNames() {
		return innerVarNames;
	}

	@Override
	public double apply(double... args) {
		return this.funCompose.apply(args);
	}
	
}
