package edu.uta.futureye.lib.shapefun;

import java.util.LinkedList;
import java.util.List;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.function.AbstractFunction;
import edu.uta.futureye.function.Variable;
import edu.uta.futureye.function.basic.FC;
import edu.uta.futureye.function.intf.Function;
import edu.uta.futureye.function.intf.ScalarShapeFunction;
import edu.uta.futureye.function.intf.ShapeFunction;
import edu.uta.futureye.util.FutureyeException;
import edu.uta.futureye.util.container.ObjList;

/**
 * 3D tetrahedral coordinates
 * 
 * 3D 四面体 局部坐标 线性形函数
 * 
 * 创建形函数需要的主要步骤：
 * 
 * 1. 创建形函数表达式（局部坐标）
 * 2. （复合）函数求值（局部坐标），（复合）函数求导数（关于物理坐标）
 * 3. 提供降维后的形函数
 * 
 * 本类采用直接定义法，并且函数值和导数已经有推导出来的表达式。
 * 复合函数法见：
 * SFQuadraticLocal2D函数值和导数没有完全推到出来，需要的计算量较大
 * SFLinearLocal2D函数值和导数已经有推导出来的表达式。
 * 
 * N = N(r,s,t,u) = N( r(x,y,z), s(x,y,z), t(x,y,z), u(x,y,z) )
 * N1 = r
 * N2 = s
 * N3 = t
 * N4 = u
 * 
 * @author liuyueming
 *
 */
public class SFLinearLocal3D extends AbstractFunction implements ScalarShapeFunction {
	private int funIndex;
	private List<String> varNames = new LinkedList<String>();
	private ObjList<String> innerVarNames = null;
	
	protected Element e = null;
	private double x1,x2,x3,x4;
	private double y1,y2,y3,y4;
	private double z1,z2,z3,z4;
	private double a1,a2,a3,a4;
	private double b1,b2,b3,b4;
	private double c1,c2,c3,c4;
	private double volume;
	
	/**
	 * 构造下列形函数中的一个：
	 * @param funID = 1,2,3,4 (N1,N2,N3,N4)
	 * 
	 */
	public void Create(int funID) {
		funIndex = funID - 1;
		if(funID<1 || funID>4) {
			System.out.println("ERROR: funID should be 1,2,3 or 4.");
			return;
		}
		
		varNames.add("r");
		varNames.add("s");
		varNames.add("t");
		varNames.add("u");
		innerVarNames = new ObjList<String>("x","y","z");
		
	}
	
	public SFLinearLocal3D(int funID) {
		this.Create(funID);
	}
	
	@Override
	public void asignElement(Element e) {
		this.e = e;
		
		x1 = e.nodes.at(1).coord(1);
		x2 = e.nodes.at(2).coord(1);
		x3 = e.nodes.at(3).coord(1);
		x4 = e.nodes.at(4).coord(1);
		y1 = e.nodes.at(1).coord(2);
		y2 = e.nodes.at(2).coord(2);
		y3 = e.nodes.at(3).coord(2);
		y4 = e.nodes.at(4).coord(2);
		z1 = e.nodes.at(1).coord(3);
		z2 = e.nodes.at(2).coord(3);
		z3 = e.nodes.at(3).coord(3);
		z4 = e.nodes.at(4).coord(3);
		
		a1=y2*(z4-z3)-y3*(z4-z2)+y4*(z3-z2);
		a2=-y1*(z4-z3)+y3*(z4-z1)-y4*(z3-z1);
		a3=y1*(z4-z2)-y2*(z4-z1)+y4*(z2-z1);
		a4=-y1*(z3-z2)+y2*(z3-z1)-y3*(z2-z1);
		
		b1=-x2*(z4-z3)+x3*(z4-z2)-x4*(z3-z2);
		b2=x1*(z4-z3)-x3*(z4-z1)+x4*(z3-z1);
		b3=-x1*(z4-z2)+x2*(z4-z1)-x4*(z2-z1);
		b4=x1*(z3-z2)-x2*(z3-z1)+x3*(z2-z1);
		
		c1=x2*(y4-y3)-x3*(y4-y2)+x4*(y3-y2);
		c2=-x1*(y4-y3)+x3*(y4-y1)-x4*(y3-y1);
		c3=x1*(y4-y2)-x2*(y4-y1)+x4*(y2-y1);
		c4=-x1*(y3-y2)+x2*(y3-y1)-x3*(y2-y1);
		
		/*
		      |x2-x1 x3-x1 x4-x1| |1 2 3|
		6*v = |y2-y1 y3-y1 y4-y1|=|4 5 6|=1*(5*9-8*6) + 4*(8*3-2*9) + 7*(2*9-8*3)
		      |z2-z1 z3-z1 z4-z1| |7 8 9|
		*/
		volume = (x2-x1)*((y3-y1)*(z4-z1)-(y4-y1)*(z3-z1))
			   + (y2-y1)*((x4-x1)*(z3-z1)-(x3-x1)*(z4-z1))
			   + (z2-z1)*((x3-x1)*(y4-y1)-(x4-x1)*(y3-y1));
		//TODO 去掉绝对值
		volume = Math.abs(volume/6.0);
		
	}

	
	SFLinearLocal2D[] faceSF = {
			new SFLinearLocal2D(1),
			new SFLinearLocal2D(2),
			new SFLinearLocal2D(3)			
		};
	@Override
	public ShapeFunction restrictTo(int funIndex) {
		return faceSF[funIndex-1];
	}

	@Override
	public Function _d(String var) {
		if(this.volume < 0.0) {
			FutureyeException e = new FutureyeException("SFLinearLocal3D: volume < 0.0");
			e.printStackTrace();
			return null;
		}
		
		//关于自由变量r,s,t求导，u为非自由变量
		if( var.equals("r") || var.equals("s") || var.equals("t") ) {
			//u=1-r-s-t
			if(funIndex == 3)
				return new FC(-1.0);
			if(var.equals(varNames.get(funIndex)))
				return new FC(1.0);
			else
				return new FC(0.0);
		} else if(var.equals("u")) {
			FutureyeException e = new FutureyeException("Error: u is not free variable");
			e.printStackTrace();
			System.exit(0);
		}
		
		//关于x,y,z求导
		if(var.equals("x")) {
			if(funIndex == 0)
				return new FC(a1/(6*volume));
			else if(funIndex == 1)
				return new FC(a2/(6*volume));
			else if(funIndex == 2)
				return new FC(a3/(6*volume));
			else if(funIndex == 3)
				return new FC(a4/(6*volume));
			else {
				FutureyeException e = new FutureyeException("Error: derivative(x)");
				e.printStackTrace();
				System.exit(0);
			}
		} else if(var.equals("y")) {
			if(funIndex == 0)
				return new FC(b1/(6*volume));
			else if(funIndex == 1)
				return new FC(b2/(6*volume));
			else if(funIndex == 2)
				return new FC(b3/(6*volume));
			else if(funIndex == 3)
				return new FC(b4/(6*volume));
			else {
				FutureyeException e = new FutureyeException("Error: derivative(y)");
				e.printStackTrace();
				System.exit(0);
			}
		} else if(var.equals("z")) {
			if(funIndex == 0)
				return new FC(c1/(6*volume));
			else if(funIndex == 1)
				return new FC(c2/(6*volume));
			else if(funIndex == 2)
				return new FC(c3/(6*volume));
			else if(funIndex == 3)
				return new FC(c4/(6*volume));
			else {
				FutureyeException e = new FutureyeException("Error: derivative(z)");
				e.printStackTrace();
				System.exit(0);
			}
		} else {
			FutureyeException e = new FutureyeException("Error: derivative()");
			e.printStackTrace();	
			System.exit(0);
		}
		return null;
	}

	@Override
	public void setVarNames(List<String> varNames) {
		this.varNames = varNames;
	}

	@Override
	public double value(Variable v) {
		if(funIndex == 0)
			return v.get("r");
		else if(funIndex == 1)
			return v.get("s");
		else if(funIndex == 2)
			return v.get("t");
		else if(funIndex == 3)
			return v.get("u");
		else {
			FutureyeException e = new FutureyeException("Error: funIndex="+funIndex);
			e.printStackTrace();
			System.exit(0);
		}
		return 0.0;
	}

	@Override
	public List<String> varNames() {
		return this.varNames;
	}
	
	public String toString() {
		return varNames.get(funIndex);
	}

	@Override
	public ObjList<String> innerVarNames() {
		return innerVarNames;
	}
}