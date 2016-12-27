package edu.uta.futureye.util;

import java.io.FileOutputStream;
import java.lang.reflect.Method;

import org.objectweb.asm.ClassWriter;
import org.objectweb.asm.Label;
import org.objectweb.asm.MethodVisitor;
import org.objectweb.asm.Opcodes;

/**
 * A wrapper class for ClassWriter This is used to generate a class
 * 
 * @author yueming.liu
 *
 */
public class ClassGenerator implements Opcodes {

	ClassWriter cw = new ClassWriter(ClassWriter.COMPUTE_MAXS
			| ClassWriter.COMPUTE_FRAMES);
	// ClassWriter cw = new ClassWriter(ClassWriter.COMPUTE_MAXS); //debug
	// purpose

	private String className;

	public ClassGenerator(String className) {
		this.className = className;
	}

	/**
	 * Start generating a class
	 * 
	 * @param interfaces
	 */
	public void startClass(String parent, String[] interfaces) {
		cw.visit(Opcodes.V1_6, ACC_PUBLIC + ACC_SUPER, this.className, null,
				parent, interfaces);

		MethodVisitor mv = cw.visitMethod(ACC_PUBLIC, "<init>", "()V", null,
				null);
		mv.visitCode();
		Label l0 = new Label();
		mv.visitLabel(l0);
		mv.visitLineNumber(3, l0);
		mv.visitVarInsn(ALOAD, 0);
		mv.visitMethodInsn(INVOKESPECIAL, parent, "<init>", "()V",
				false);
		mv.visitInsn(RETURN);
		Label l1 = new Label();
		mv.visitLabel(l1);
		mv.visitLocalVariable("this", "L" + this.className + ";", null, l0, l1,
				0);
		mv.visitMaxs(1, 1);
		mv.visitEnd();
	}

	/**
	 * End generating a class
	 */
	public void endClass() {
		cw.visitEnd();
	}

	/**
	 * Start generating a method of the class
	 * 
	 * @param access
	 * @param name
	 * @param type
	 * @return
	 */
	public MethodVisitor startMethod(int access, String name, String type) {
		return cw.visitMethod(access, name, type, null, null);
	}

	/**
	 * Start generating code of a method
	 * 
	 * @param mv
	 * @param startLabel
	 */
	public void startCode(MethodVisitor mv, Label startLabel) {
		mv.visitCode();
		mv.visitLabel(startLabel);
	}

	/**
	 * Get ClassWriter
	 * 
	 * @return
	 */
	public ClassWriter getClassWriter() {
		return this.cw;
	}

	/**
	 * End generating code of a method
	 * 
	 * @param mv
	 * @param endLabel
	 */
	public void endCode(MethodVisitor mv, Label endLabel) {
		mv.visitLabel(endLabel);
		mv.visitEnd();
	}

	/**
	 * Dump the byte array of the generated class
	 * 
	 * @return
	 * @throws Exception
	 */
	public byte[] dump() throws Exception {
		return cw.toByteArray();
	}

	public String getClassName() {
		return this.className;
	}

	/**
	 * A simple example use `javap -v Test.class` to print the generated human
	 * readable bytecode
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		FuncClassLoader mcl = new FuncClassLoader(
				ClassGenerator.class.getClassLoader());
		try {
			ClassGenerator cgen = new ClassGenerator("Test");
			cgen.startClass("java/lang/Object", null);

			Label startLabel = new Label();
			Label endLabel = new Label();
			MethodVisitor mv = cgen.startMethod(Opcodes.ACC_PUBLIC, "eval",
					"(DD)D");
			cgen.startCode(mv, startLabel);

			// Generate 10.0+2.0
			mv.visitLdcInsn(10.0);
			mv.visitLdcInsn(2.0);
			mv.visitInsn(DADD);
			mv.visitInsn(DRETURN);

			mv.visitLocalVariable("this", "L" + cgen.getClassName() + ";",
					null, startLabel, endLabel, 0);
			mv.visitLocalVariable("a", "D", null, startLabel, endLabel, 1);
			mv.visitLocalVariable("b", "D", null, startLabel, endLabel, 3);
			mv.visitMaxs(4, 5);

			cgen.endCode(mv, endLabel);
			cgen.endClass();

			byte[] bcode = cgen.dump();
			Class<?> c = mcl.defineClassForName(null, bcode);
			Method m1 = c.getMethod("eval", double.class, double.class);
			Object o = c.newInstance();
			// The parameter is unused in the generated class
			System.out.println(m1.invoke(o, 1, 2));

			FileOutputStream fos = new FileOutputStream(cgen.getClassName()
					+ ".class");
			fos.write(bcode);
			fos.close();

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static String getASMName(Class<?> c) {
		return c.getName().replaceAll("\\.", "/");
	}
}
