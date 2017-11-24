/*

  TEST Java Generator

  create:
   includes files .h
   sources .cpp
   make build.xml

*/

package com.numx.automatization.generator;

import java.io.File;
import java.io.FileWriter;
import java.util.Calendar;
import java.util.GregorianCalendar;
import java.util.Iterator;
import java.util.List;

import org.apache.velocity.VelocityContext;
import org.apache.velocity.context.Context;

import com.numx.automatization.model.Function;
import com.numx.automatization.model.Module;
import com.numx.automatization.model.Product;
import com.numx.automatization.xml.XMLParser;

public class Java_TEST_Generator extends VelocityHelper {

	public final static String TEMPLATE_TEST_JAVA_DIR = "com/numx/automatization/generator/test";
	
	private Context ctx;
	
	/**
	 * @param args
	 * args[0] is the path of descriptor files
	 * args[1] is the path of root of NumX
	 */
	public static void main(String[] args) throws Exception 
	{
		//reading descripor files and getting a list of elements of type Product
		XMLParser parser = new XMLParser();
		List products = parser.getProducts(args[0]);
		
		//generating files
		String outDir = args[1];
		System.out.println("Generating tests for Java ...");
		new Java_TEST_Generator().generate(products, outDir);
	}

	public Java_TEST_Generator() throws Exception{
		super();
		VMs.add(TEMPLATE_TEST_JAVA_DIR);
        initVelocityEngine();
	}

	public void generate(List products, String outDir) throws Exception {
        ctx = new VelocityContext();
        ctx.put("year", new Integer((new GregorianCalendar()).get(Calendar.YEAR)).toString());

   	for (Iterator itProds = products.iterator(); itProds.hasNext();) {
			Product prod = (Product) itProds.next();
	        ctx.put("prod", prod);
	        
	        String outTestDir = outDir + "/test";
	        generateTestByProduct(prod, outTestDir);

			int index = 10;
        	for (Iterator itMods = prod.getModules().iterator(); itMods.hasNext();) {
        		Module mod = (Module) itMods.next();
        		if (!mod.is_onlyC() && mod.is_release()) {        			
	    	        ctx.put("mod", mod);
					ctx.put("LANG", "java");
					ctx.put("ModProjectGUIDJAVA", "9F113DAA-866E-4870-B356-EB23C0B57A"+Integer.toString(index++));

					generateTestByModule(prod, mod, outTestDir);
	
	            	for (Iterator itFuncs = mod.getFunctions().iterator(); itFuncs.hasNext();) {
	            		Function func = (Function) itFuncs.next();
	            		if (!func.is_onlyC() && func.is_release()) {        		
		        	        ctx.put("func", func);

							if (!("NumX" + mod.get_moduleNameLowerNoSpace() + "Test").equals(func.get_functionName())) {
		        				generateTestByFunction(func, outTestDir);
							}
	            		}
	            	}
        		}
        	}
    	}
	}

	public void generateTestByProduct(Product prod, String outTestDir) throws Exception {

        File f = new File(outTestDir + "/src/java", "Test" + prod.get_productShortName() + ".java");
		f.getParentFile().mkdirs();
        FileWriter fw = new FileWriter(f);
        velocityEngine.getTemplate(TEMPLATE_TEST_JAVA_DIR +"/ProductJava.vm").merge(ctx, fw);
        fw.flush();
        fw.close();	        

        f = new File(outTestDir, "build.xml");
		f.getParentFile().mkdirs();
        fw = new FileWriter(f);
        velocityEngine.getTemplate(TEMPLATE_TEST_JAVA_DIR +"/Build.vm").merge(ctx, fw);
        fw.flush();
        fw.close();	        
	}
	
	public void generateTestByModule(Product prod, Module mod, String outTestDir) throws Exception {

		File f = new File(outTestDir + "/src/java", "Test" + prod.get_productShortName() + mod.get_moduleNameNoSpace() + ".java");
		f.getParentFile().mkdirs();
        FileWriter fw = new FileWriter(f);
        velocityEngine.getTemplate(TEMPLATE_TEST_JAVA_DIR +"/ModuleJava.vm").merge(ctx, fw);
        fw.flush();
        fw.close();

	}

	public void generateTestByFunction(Function func, String outTestDir) throws Exception {

		File f = new File(outTestDir + "/src/java/", "Test" + func.get_functionName() + ".java");
		f.getParentFile().mkdirs();
        FileWriter fw = new FileWriter(f);
        velocityEngine.getTemplate(TEMPLATE_TEST_JAVA_DIR +"/FunctionJava.vm").merge(ctx, fw);
        fw.flush();
        fw.close();

	}
	
}
