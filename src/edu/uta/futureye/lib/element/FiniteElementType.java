package edu.uta.futureye.lib.element;

import edu.uta.futureye.core.Element;
import edu.uta.futureye.core.Mesh;

public interface FiniteElementType {
	/**
	 * Associate degrees of freedom (DOF) to element e
	 * 
	 * @param e
	 */
	void assignTo(Element e);
	
	/**
	 * 初始化自由度编号生成器
	 * @param nTotalNodes
	 */
	void initDOFIndexGenerator(Mesh mesh);
	
	/**
	 * 获得向量值形函数的维度
	 * @return
	 */
	int getVectorShapeFunctionDim();
	
	/**
	 * 获得单元上，自由度总数。如果单元上自由度关联的形函数为标量函数，
	 * 参数vsfDim可取任意值，如果为向量值函数，需要指定向量维度vsfDim，
	 * 将返回该维度对应的自由度总数。
	 * 
	 * @param vsfDim 向量值形函数的维度
	 * @return
	 */
	int getDOFNumOnElement(int vsfDim);
	
	/**
	 * 获得整个网格上，自由度总数。如果网格的单元上自由度关联的形函数为标量函数，
	 * 参数vsfDim可取任意值，如果为向量值函数，需要指定向量维度vsfDim，
	 * 将返回该维度对应的自由度总数。
	 * 
	 * @param mesh
	 * @param vsfDim
	 * @return
	 */
	int getDOFNumOnMesh(Mesh mesh,int vsfDim);
}
