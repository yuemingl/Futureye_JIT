package edu.uta.futureye.function.basic;

import java.util.ArrayList;
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
import edu.uta.futureye.function.AbstractMathFunc;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.VariableArray;
import edu.uta.futureye.function.intf.MathFunc;
import edu.uta.futureye.util.BytecodeUtils;
import edu.uta.futureye.util.FuncClassLoader;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.Utils;

/**
 * Composite function
 * <p><blockquote><pre>
 * For example:
 *  MathFunc f = r*s + 1;
 *  Map<String, MathFunc> fInners = new HashMap<String, MathFunc>();
 *  fInners.put("r", x*x);
 *  fInners.put("s", y+1);
 *  MathFunc fc = f.compose(fInners);
 *  System.out.println(fc); //f(x,y) = (x*x)*(y + 1.0) + 1.0
 * </pre></blockquote>
 */
public class FComposite extends AbstractMathFunc {
	public MathFunc fOuter;
	public Map<String,MathFunc> fInners;
	boolean isOuterVariablesActive;
	
	public FComposite(MathFunc fOuter, Map<String,MathFunc> fInners) {
		this.fOuter = fOuter;
		
		//Extends variable names in fInners (copy on change)
		List<String> list = new ArrayList<String>();
		for(Entry<String, MathFunc> e : fInners.entrySet()) {
			list = Utils.mergeList(list, e.getValue().getVarNames());
		}
		Map<String,MathFunc> fInners2 = new HashMap<String,MathFunc>();
		Map<String, Integer> argsMap = Utils.getIndexMap(list);
		for(Entry<String, MathFunc> e : fInners.entrySet()) {
			if(!Utils.isMapContain(argsMap, e.getValue().getArgIdxMap())) {
				MathFunc f = e.getValue().copy().setArgIdx(argsMap);
				fInners2.put(e.getKey(), f);
			} else {
				fInners2.put(e.getKey(), e.getValue());
			}
		}
		this.fInners = fInners2;
		
		// Default to use free variables in fInners
		this.setVarNames(list);
		this.setArgIdx(Utils.getIndexMap(list));
		isOuterVariablesActive = false;
		//this.setVarNames(fOuter.getVarNames());
		//this.setArgIdx(Utils.getIndexMap(fOuter.getVarNames()));
	}

	@Override
	public double apply(Variable v) {
		return apply(v,null);
	}

	@Override
	public double apply(double... args) {
		return apply(null, null, args);
	}

	@Override
	public double apply(Element e, Node n, double... args) {
		if(this.isOuterVariablesActive) {
			return fOuter.apply(e, n, args);
		} else {
			List<String> vn = fOuter.getVarNames();
			double[] newArgs = new double[vn.size()];
			for(int i=0; i<vn.size(); i++) {
				MathFunc fInner = fInners.get(vn.get(i));
				newArgs[i] = fInner.apply(e, n, args);
			}
			return fOuter.apply(newArgs);
		}
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
	public MathFunc diff(String varName) {
		MathFunc rlt = null;
		if(fOuter.getVarNames().contains(varName)) {
			//f(x,y)关于x或y求导
			rlt = fOuter.diff(varName);
			return rlt;
		} else {
			//f(x,y)关于r或s求导
			rlt = new FC(0.0);
			for(String innerVarName : fOuter.getVarNames()) {
				MathFunc fInner = fInners.get(innerVarName);
				if(fInner != null) {
					MathFunc rltOuter = fOuter.diff(innerVarName);
					if(!(rltOuter.isConstant()))
						rltOuter = rltOuter.compose(fInners);
					MathFunc rltInner = fInner.diff(varName);
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
	public String getExpr() {
		if(this.isOuterVariablesActive) {
			return fOuter.getExpr();
		} else {
			String rlt = fOuter.getExpr();
			for(Entry<String,MathFunc> map : fInners.entrySet()) {
	//			String names = map.getValue().getVarNames().toString();
	//			rlt = rlt.replace(map.getKey(), 
	//					map.getKey()+"("+names.substring(1,names.length()-1)+")");
				if(map.getValue().getOpOrder() == OP_ORDER0)
					rlt = rlt.replace(map.getKey(), map.getValue().getExpr());
				else
					rlt = rlt.replace(map.getKey(), "("+map.getValue().getExpr()+")");
			}
			return rlt;
		}
	}

	@Override
	public InstructionHandle bytecodeGen(String clsName, MethodGen mg,
			ConstantPoolGen cp, InstructionFactory factory,
			InstructionList il, Map<String, Integer> argsMap, int argsStartPos, 
			Map<MathFunc, Integer> funcRefsMap) {
		if(this.isOuterVariablesActive) {
			return fOuter.bytecodeGen(clsName, mg, cp, factory, il, argsMap, argsStartPos, funcRefsMap);
		} else {
			// Prepare a double array as the arguments for the fOuter function
			// which is equal to call the outer function
			LocalVariableGen lg;
			
			//double[] arg = null;
			lg = mg.addLocalVariable("aryArgOuter",
				new ArrayType(Type.DOUBLE, 1), null, null);
			int aryArgOuter = lg.getIndex();
			il.append(InstructionConstants.ACONST_NULL);
			lg.setStart(il.append(new ASTORE(aryArgOuter))); // "idxArg" valid from here
			
			//arg = new double[size]
			//il.append(new PUSH(cp, fInners.size()));
			il.append(new PUSH(cp, fOuter.getVarNames().size()));
			il.append(new NEWARRAY(Type.DOUBLE));
			il.append(new ASTORE(aryArgOuter));
			
			//int index = 0;
			Map<String, Integer> argMap = fOuter.getArgIdxMap();
			for(String name : fOuter.getVarNames()) {
				MathFunc f = fInners.get(name);
				//aryArgOuter[argIdx] = {value of the argument}
				il.append(new ALOAD(aryArgOuter));
				il.append(new PUSH(cp, argMap.get(name))); //index++
				if(f != null) {
					List<String> args = f.getVarNames();
					HashMap<String, Integer> fArgsMap = new HashMap<String, Integer>();
					for(int i=0; i<args.size(); i++) {
						fArgsMap.put(args[i], argsMap.get(args[i]));
					}
					f.bytecodeGen(clsName, mg, cp, factory, il, fArgsMap, 3, funcRefsMap);
				} else {
					il.append(new PUSH(cp, 0.0)); //pad 0.0
				}
				il.append(new DASTORE());
			}
			// Pass a double array to fOuter
			return fOuter.bytecodeGen(clsName, mg, cp, factory, il, fOuter.getArgIdxMap(), aryArgOuter, funcRefsMap);
		}
	}
	
	public InstructionHandle bytecodeGenSlow(String clsName, MethodGen mg,
			ConstantPoolGen cp, InstructionFactory factory,
			InstructionList il, Map<String, Integer> argsMap, int argsStartPos, 
			Map<MathFunc, Integer> funcRefsMap) {
		if(this.isOuterVariablesActive) {
			return fOuter.bytecodeGen(clsName, mg, cp, factory, il, argsMap, argsStartPos, funcRefsMap);
		} else {
			String outerName  = "fun_outer_"+java.util.UUID.randomUUID().toString().replaceAll("-", "");
			// Generate the outer function
			FuncClassLoader<CompiledFunc> fcl = new FuncClassLoader<CompiledFunc>();
			ClassGen genClass = BytecodeUtils.genClass(fOuter, null, outerName, true, true);
			fcl.newInstance(genClass);
	
			// Prepare arguments for calling the outer function
			LocalVariableGen lg;
			//double[] arg = null;
			lg = mg.addLocalVariable("arg_"+outerName,
				new ArrayType(Type.DOUBLE, 1), null, null);
			int idxArg = lg.getIndex();
			il.append(InstructionConstants.ACONST_NULL);
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
					fArgsMap.put(args[i], argsMap.get(args[i]));
				}
				f.bytecodeGen(clsName, mg, cp, factory, il, fArgsMap, 3, funcRefsMap);
				il.append(new DASTORE());
			}
			
			// Call the outer function
			il.append(InstructionConstants.ACONST_NULL);
			il.append(InstructionConstants.ACONST_NULL);
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
	
	@Override
	public MathFunc setArgIdx(Map<String, Integer> argsMap) {
		super.setArgIdx(argsMap);
		if(this.isOuterVariablesActive) {
			this.fOuter.setArgIdx(argsMap);
		} else {
		}
		return this;
	}
	
	@Override
	public MathFunc setActiveVarNames(List<String> varNames) {
		if(Utils.isListEqualIgnoreOrder(fOuter.getVarNames(),  varNames)) {
			this.isOuterVariablesActive = true;
			this.setVarNames(varNames);
			this.setArgIdx(Utils.getIndexMap(varNames));
			return this;
		} else {
			List<String> list = new ArrayList<String>();
			for(Entry<String, MathFunc> e : fInners.entrySet()) {
				list = Utils.mergeList(list, e.getValue().getVarNames());
			}
			if(Utils.isListEqualIgnoreOrder(list, varNames)){
				this.isOuterVariablesActive = false;
				this.setVarNames(varNames);
				Map<String, Integer> argsMap = Utils.getIndexMap(list);
				for(Entry<String, MathFunc> e : fInners.entrySet()) {
					e.getValue().setArgIdx(argsMap);
				}
				this.setArgIdx(Utils.getIndexMap(list));
				return this;
			}
		}
		throw new FutureyeException("Active variable names are different from all existing variable names!");
	}
}
