package edu.uta.futureye.test.junit;

import static org.junit.Assert.*;

import org.junit.Test;

import edu.uta.futureye.bytecode.CompiledFunc;
import edu.uta.futureye.function.basic.FXY;

public class TestFunctions {

	@Test
	public void testFXY() {
		FXY f = new FXY(2,3,4); //2x+3y+4
		CompiledFunc cf = f.compileWithASM(new String[]{"x","y"});
		double d = cf.apply(new double[]{3,3});
		assertTrue(Math.abs(d-19)<1e-8);
	}


}
