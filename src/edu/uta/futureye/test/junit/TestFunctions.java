package edu.uta.futureye.test.junit;

import static org.junit.Assert.*;
import static edu.uta.futureye.function.FMath.*;

import org.junit.Test;

import edu.uta.futureye.bytecode.CompiledFunc;
import edu.uta.futureye.function.FMath;
import edu.uta.futureye.function.basic.FXY;

public class TestFunctions {
	@Test
	public void testFXY() {
		FXY f = new FXY(2,3,4); //2x+3y+4
		CompiledFunc cf = f.compileWithASM(new String[]{"x","y"});
		double d = cf.apply(new double[]{3,3});
		assertTrue(Math.abs(d-19)<1e-8);
		
		assertTrue(Math.abs(Math.sin(1.0)     -FMath.sin(x)  .compileWithASM(new String[]{"x"}).apply(1.0))<1e-8);
		assertTrue(Math.abs(Math.cos(1.0)     -FMath.cos(x)  .compileWithASM(new String[]{"x"}).apply(1.0))<1e-8);
		assertTrue(Math.abs(Math.asin(1.0)    -FMath.asin(x) .compileWithASM(new String[]{"x"}).apply(1.0))<1e-8);
		assertTrue(Math.abs(Math.acos(1.0)    -FMath.acos(x) .compileWithASM(new String[]{"x"}).apply(1.0))<1e-8);
		assertTrue(Math.abs(Math.sinh(1.0)    -FMath.sinh(x) .compileWithASM(new String[]{"x"}).apply(1.0))<1e-8);
		assertTrue(Math.abs(Math.cosh(1.0)    -FMath.cosh(x) .compileWithASM(new String[]{"x"}).apply(1.0))<1e-8);
		assertTrue(Math.abs(Math.tan(1.0)     -FMath.tan(x)  .compileWithASM(new String[]{"x"}).apply(1.0))<1e-8);
		assertTrue(Math.abs(Math.tanh(1.0)    -FMath.tanh(x) .compileWithASM(new String[]{"x"}).apply(1.0))<1e-8);
		assertTrue(Math.abs(Math.log(1.0)     -FMath.log(x)  .compileWithASM(new String[]{"x"}).apply(1.0))<1e-8);
		assertTrue(Math.abs(Math.exp(2.0)     -FMath.exp(x)  .compileWithASM(new String[]{"x"}).apply(2.0))<1e-8);
		assertTrue(Math.abs(Math.log10(1.0)   -FMath.log10(x).compileWithASM(new String[]{"x"}).apply(1.0))<1e-8);
		assertTrue(Math.abs(Math.pow(2.0, 5)  -FMath.pow(x,y).compileWithASM(new String[]{"x","y"}).apply(2.0, 5))<1e-8);
		assertTrue(Math.abs(Math.pow(2.0, 3.3)-FMath.pow(x,y).compileWithASM(new String[]{"x","y"}).apply(2.0, 3.3))<1e-8);
		assertTrue(Math.abs(Math.signum(1.0)  -FMath.signum(x).compileWithASM(new String[]{"x"}).apply(1.0))<1e-8);
		assertTrue(Math.abs(Math.abs(-1.0)     -FMath.abs(x)  .compileWithASM(new String[]{"x"}).apply(-1.0))<1e-8);
		assertTrue(Math.abs(Math.sqrt(1.0)    -FMath.sqrt(x) .compileWithASM(new String[]{"x"}).apply(1.0))<1e-8);
		assertTrue(Math.abs(Math.max(3.0, 5.0)-FMath.max(x,y).compileWithASM(new String[]{"x","y"}).apply(3.0, 5.0))<1e-8);
		assertTrue(Math.abs(Math.min(3.0, 5.0)-FMath.min(x,y).compileWithASM(new String[]{"x","y"}).apply(3.0, 5.0))<1e-8);
	}


}
