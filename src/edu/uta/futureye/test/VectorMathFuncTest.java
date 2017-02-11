package edu.uta.futureye.test;

import org.junit.Test;

import edu.uta.futureye.function.basic.SpaceVectorFunction;
import edu.uta.futureye.function.intf.VectorMathFunc;
import static edu.uta.futureye.function.FMath.*;
import static org.junit.Assert.*;

public class VectorMathFuncTest {
	static double eps = 1e-6;
	@Test
	public void testBasic() {
		VectorMathFunc f = new SpaceVectorFunction(x,y,z);
		VectorMathFunc g = new SpaceVectorFunction(C(2),C(2),C(2));
		
		assertTrue(Math.abs(12.0-dot(f,g).apply(1,2,3))<eps);

		VectorMathFunc h = null;
		
		h = f + g;
		assertTrue(Math.abs(12.0-h[1].apply(10))<eps);
		assertTrue(Math.abs(12.0-h[2].apply(10))<eps);
		assertTrue(Math.abs(12.0-h[3].apply(10))<eps);
		
		h = f - g;
		assertTrue(Math.abs(8.0-h[1].apply(10))<eps);
		assertTrue(Math.abs(8.0-h[2].apply(10))<eps);
		assertTrue(Math.abs(8.0-h[3].apply(10))<eps);
		
		h = f * g;
		assertTrue(Math.abs(20.0-h[1].apply(10))<eps);
		assertTrue(Math.abs(20.0-h[2].apply(10))<eps);
		assertTrue(Math.abs(20.0-h[3].apply(10))<eps);

		h = f / g;
		assertTrue(Math.abs(5.0-h[1].apply(10))<eps);
		assertTrue(Math.abs(5.0-h[2].apply(10))<eps);
		assertTrue(Math.abs(5.0-h[3].apply(10))<eps);

		h = -f;
		assertTrue(Math.abs(-10.0-h[1].apply(10))<eps);
		assertTrue(Math.abs(-10.0-h[2].apply(10))<eps);
		assertTrue(Math.abs(-10.0-h[3].apply(10))<eps);
		
		h = f + 2;
		assertTrue(Math.abs(12.0-h[1].apply(10))<eps);
		assertTrue(Math.abs(12.0-h[2].apply(10))<eps);
		assertTrue(Math.abs(12.0-h[3].apply(10))<eps);
		h = 2 + f;
		assertTrue(Math.abs(12.0-h[1].apply(10))<eps);
		assertTrue(Math.abs(12.0-h[2].apply(10))<eps);
		assertTrue(Math.abs(12.0-h[3].apply(10))<eps);

		h = f - 2;
		assertTrue(Math.abs(8.0-h[1].apply(10))<eps);
		assertTrue(Math.abs(8.0-h[2].apply(10))<eps);
		assertTrue(Math.abs(8.0-h[3].apply(10))<eps);
		h = 2 - f;
		assertTrue(Math.abs(-8.0-h[1].apply(10))<eps);
		assertTrue(Math.abs(-8.0-h[2].apply(10))<eps);
		assertTrue(Math.abs(-8.0-h[3].apply(10))<eps);

		h = f * 2;
		assertTrue(Math.abs(20.0-h[1].apply(10))<eps);
		assertTrue(Math.abs(20.0-h[2].apply(10))<eps);
		assertTrue(Math.abs(20.0-h[3].apply(10))<eps);

		h = 2 * f;
		assertTrue(Math.abs(20.0-h[1].apply(10))<eps);
		assertTrue(Math.abs(20.0-h[2].apply(10))<eps);
		assertTrue(Math.abs(20.0-h[3].apply(10))<eps);

		h = f / 2;
		assertTrue(Math.abs(5.0-h[1].apply(10))<eps);
		assertTrue(Math.abs(5.0-h[2].apply(10))<eps);
		assertTrue(Math.abs(5.0-h[3].apply(10))<eps);
		
		h = 2 / f;
		assertTrue(Math.abs(2.0/10.0-h[1].apply(10))<eps);
		assertTrue(Math.abs(2.0/10.0-h[2].apply(10))<eps);
		assertTrue(Math.abs(2.0/10.0-h[3].apply(10))<eps);
	}
}
