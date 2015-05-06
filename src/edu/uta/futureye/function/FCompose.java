package edu.uta.futureye.function;

import static com.sun.org.apache.bcel.internal.generic.InstructionConstants.ACONST_NULL;
import static com.sun.org.apache.bcel.internal.generic.InstructionConstants.DRETURN;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import com.sun.org.apache.bcel.internal.Constants;
import com.sun.org.apache.bcel.internal.generic.ALOAD;
import com.sun.org.apache.bcel.internal.generic.ASTORE;
import com.sun.org.apache.bcel.internal.generic.ArrayType;
import com.sun.org.apache.bcel.internal.generic.ClassGen;
import com.sun.org.apache.bcel.internal.generic.ConstantPoolGen;
import com.sun.org.apache.bcel.internal.generic.DASTORE;
import com.sun.org.apache.bcel.internal.generic.InstructionConstants;
import com.sun.org.apache.bcel.internal.generic.InstructionFactory;
import com.sun.org.apache.bcel.internal.generic.InstructionHandle;
import com.sun.org.apache.bcel.internal.generic.InstructionList;
import com.sun.org.apache.bcel.internal.generic.LocalVariableGen;
import com.sun.org.apache.bcel.internal.generic.MethodGen;
import com.sun.org.apache.bcel.internal.generic.NEWARRAY;
import com.sun.org.apache.bcel.internal.generic.PUSH;
import com.sun.org.apache.bcel.internal.generic.Type;

import edu.uta.futureye.bytecode.CompiledFunc;
import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Node;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.util.BytecodeUtils;
import edu.uta.futureye.util.FuncClassLoader;

public class FCompose extends AbstractMathFun {
	public MathFunc fOuter;
	public Map<String,MathFunc> fInners;
	
	public FCompose(MathFunc fOuter, Map<String,MathFunc> fInners) {
		this.fOuter = fOuter;
		this.fInners = fInners;
		this.setVarNames(fOuter.getVarNames());
	}

	@Override
	public double apply(Variable v) {
		return apply(v,null);
	}

	@Override
	public double apply(Element e, Node n, double... args) {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public double apply(Variable v, Map<Object,Object> cache) {
		
		//if(fOuter.varNames().size() == 0) {
		//	throw new FutureyeException("\nERROR:\n fOuter varNames list is empty!");
		//}
		
		//bugfix 增加或条件 
		//bugfix 3/19/12
		//bug?3/20/12  v=[r], fOuter.varNames()=[s,t], 但fOuter的表达式只有r, 这种情况下会进入else分支，
		//一般来说是不会有这种情况的，如果确实有这种情况，需要在函数类增加activeVarNames
		//if(fOuter.varNames().containsAll(v.getValues().keySet()) ||
		//		v.getValues().keySet().containsAll(fOuter.varNames())) {
		if(v.getNameValuePairs().keySet().containsAll(fOuter.getVarNames())) {
			return fOuter.apply(v,cache);
		//} else if(fOuter.varNames().size() == fInners.size()){
		} else {
			Variable newVar = new Variable();
			for(String varName : fOuter.getVarNames()) {
				MathFunc fInner = fInners.get(varName);
				if(fInner != null ) 
					newVar.set(varName, fInner.apply(v,cache));
				else //for mixed case: fOuter( x(r,s,t), y(r,s,t), r, s) bugfix 3/19/12
					newVar.set(varName, v.get(varName));
				//	throw new FutureyeException("\nERROR:\n Can not find "+varName+" in fInners.");
			}
			return fOuter.apply(newVar,cache);
		}
//		else {
//			throw new FutureyeException(
//					"\nERROR:\n Variable number mismatch of fOuter("+
//					fOuter.varNames()+") and fInner("+fInners+").");
//		}
	} 

	@Override
	public double[] applyAll(VariableArray v, Map<Object,Object> cache) {
		//bugfix 增加或条件
		if(v.getValues().keySet().containsAll(fOuter.getVarNames())) {
			return fOuter.applyAll(v,cache);
		} else {
			VariableArray newVar = new VariableArray();
			for(String varName : fOuter.getVarNames()) {
				MathFunc fInner = fInners.get(varName);
				if(fInner != null )
					newVar.set(varName, fInner.applyAll(v,cache));
				else //for mixed case: fOuter( x(r,s,t), y(r,s,t), r, s)
					newVar.set(varName, v.get(varName));
			}
			return fOuter.applyAll(newVar,cache);
		}
	}
	
	/**
	 * 链式求导
	 * f( x(r,s),y(r,s) )_r = f_x * x_r + f_y * y_r
	 */
	@Override
	public MathFunc _d(String varName) {
		MathFunc rlt = null;
		if(fOuter.getVarNames().contains(varName)) {
			//f(x,y)关于x或y求导
			rlt = fOuter._d(varName);
			return rlt;
		} else {
			//f(x,y)关于r或s求导
			rlt = new FC(0.0);
			for(String innerVarName : fOuter.getVarNames()) {
				MathFunc fInner = fInners.get(innerVarName);
				if(fInner != null) {
					MathFunc rltOuter = fOuter._d(innerVarName);
					if(!(rltOuter.isConstant()))
						rltOuter = rltOuter.compose(fInners);
					MathFunc rltInner = fInner._d(varName);
					//f_x * x_r + f_y * y_r
					rlt = rlt.A(
							rltOuter.M(rltInner)
							);
				}
			}
			return rlt;
		}
	}
	
	@Override
	public int getOpOrder() {
		return fOuter.getOpOrder();
	}
	
	@Override
	public String toString() {
		String rlt = fOuter.toString();
		for(Entry<String,MathFunc> map : fInners.entrySet()) {
			String names = map.getValue().getVarNames().toString();
			rlt = rlt.replace(map.getKey(), 
					map.getKey()+"("+names.substring(1,names.length()-1)+")");
		}
		return rlt;
	}

	@Override
	public InstructionHandle bytecodeGen(String clsName, MethodGen mg,
			ConstantPoolGen cp, InstructionFactory factory,
			InstructionList il, Map<String, Integer> argsMap, int argsStartPos, Map<MathFunc, Integer> funcRefsMap) {
		String outerName  = "fun_outer_"+java.util.UUID.randomUUID().toString().replaceAll("-", "");
		// Generate the outer function
		FuncClassLoader<CompiledFunc> fcl = new FuncClassLoader<CompiledFunc>();
		ClassGen genClass = BytecodeUtils.genClass(fOuter, outerName, true, true);
		fcl.newInstance(genClass);

		// Prepare arguments for calling the outer function
		LocalVariableGen lg;
		//double[] arg = null;
		lg = mg.addLocalVariable("arg_"+outerName,
			new ArrayType(Type.DOUBLE, 1), null, null);
		int idxArg = lg.getIndex();
		il.append(ACONST_NULL);
		lg.setStart(il.append(new ASTORE(idxArg))); // "idxArg" valid from here
		//arg = new double[size]
		il.append(new PUSH(cp, fInners.size()));
		il.append(new NEWARRAY(Type.DOUBLE));
		il.append(new ASTORE(idxArg));
		
		int index = 0;
		for(String name : fOuter.getVarNames()) {
			il.append(new ALOAD(idxArg));
			il.append(new PUSH(cp, index++));
			MathFunc f = fInners.get(name);
			HashMap<String, Integer> fArgsMap = new HashMap<String, Integer>();
			List<String> args = f.getVarNames();
			for(int i=0; i<args.size(); i++) {
				fArgsMap.put(args[i], i);
			}
			f.bytecodeGen(null, mg, cp, factory, il, fArgsMap, 3, null);
			il.append(new DASTORE());
		}
		
		// Call the outer function
		il.append(ACONST_NULL);
		il.append(ACONST_NULL);
		il.append(new ALOAD(idxArg));
		return  il.append(factory.createInvoke("edu.uta.futureye.bytecode."+outerName, "apply",
				Type.DOUBLE, 
				new Type[] { 
					Type.getType(Element.class),
					Type.getType(Node.class),
					new ArrayType(Type.DOUBLE, 1)
				}, 
		Constants.INVOKESTATIC));
	}

}
