/*

 TEST C Generator

 create:
  includes files .h
  sources .cpp

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

public class C_TEST_Generator extends VelocityHelper {

	// paths of velocity models
	public final static String TEMPLATE_TEST_CPP_DIR = "com/numx/automatization/generator/test";
	
	private Context ctx;
	
	/**
	 * @param args
	 * args[0] is the path of descriptor files
	 * args[1] is the path of root of NumX
	 */
	public static void main(String[] args) throws Exception{
		//reading descripor files and getting a list of elements of type Product
		XMLParser parser = new XMLParser();
		List products = parser.getProducts(args[0]);

		//generating files
		String outDir = args[1];
		System.out.println("Generating tests for C/C++ ...");
		new C_TEST_Generator().generate(products, outDir);
	}

	public C_TEST_Generator() throws Exception{
		super();
		VMs.add(TEMPLATE_TEST_CPP_DIR);
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
        		if (mod.is_release()) {
	    	        ctx.put("mod", mod);
					ctx.put("LANG", "c");
					ctx.put("ModProjectGUIDC", "6709D3AC-ED44-43E4-BE68-4DB3D8F4DC"+Integer.toString(index++));
					
	    	        generateTestByModule(prod, mod, outTestDir);
	
	            	for (Iterator itFuncs = mod.getFunctions().iterator(); itFuncs.hasNext();) {
	            		Function func = (Function) itFuncs.next();
	            		if (func.is_release()) {        		
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
	
	// TESTS C --- creation des fichiers product pour les tests
	public void generateTestByProduct(Product prod, String outTestDir) throws Exception {

		// fichiers tests .cpp
        File f = new File(outTestDir + "/src/cpp", "Test" + prod.get_productShortName() + ".cpp");
		f.getParentFile().mkdirs();
        FileWriter fw = new FileWriter(f);
        velocityEngine.getTemplate(TEMPLATE_TEST_CPP_DIR +"/ProductC.vm").merge(ctx, fw);
        fw.flush();
        fw.close();	        

		// ajout du fichier stdafx.h
        f = new File(outTestDir + "/src/cpp", "stdafx.h");
		f.getParentFile().mkdirs();
        fw = new FileWriter(f);
        velocityEngine.getTemplate(TEMPLATE_TEST_CPP_DIR +"/Stdafx.vm").merge(ctx, fw);
        fw.flush();
        fw.close(); 

	}

	// TESTS C -- creation des fichiers modules .cpp .h pour les tests 
	public void generateTestByModule(Product prod, Module mod, String outTestDir) throws Exception {

		File f = new File(outTestDir + "/src/cpp", "Test" + prod.get_productShortName() + mod.get_moduleNameNoSpace() + ".cpp");
		f.getParentFile().mkdirs();
        FileWriter fw = new FileWriter(f);
        velocityEngine.getTemplate(TEMPLATE_TEST_CPP_DIR +"/ModuleC.vm").merge(ctx, fw);
        fw.flush();
        fw.close();

		f = new File(outTestDir + "/src/cpp", "Test" + prod.get_productShortName() + mod.get_moduleNameNoSpace() + ".h");
		f.getParentFile().mkdirs();
        fw = new FileWriter(f);
        velocityEngine.getTemplate(TEMPLATE_TEST_CPP_DIR +"/ModuleH.vm").merge(ctx, fw);
        fw.flush();
        fw.close();

	}

	// TESTS C -- creation des fichiers fonctions .cpp .h pour les tests
	public void generateTestByFunction(Function func, String outTestDir) throws Exception {

		File f = new File(outTestDir + "/src/cpp", "Test" + func.get_functionName() + ".cpp");
		f.getParentFile().mkdirs();
        FileWriter fw = new FileWriter(f);
        velocityEngine.getTemplate(TEMPLATE_TEST_CPP_DIR +"/FunctionC.vm").merge(ctx, fw);
        fw.flush();
        fw.close();

		//System.out.printf("\n"+func.get_functionName());
		f = new File(outTestDir + "/src/cpp", "Test" + func.get_functionName() + ".h");
		f.getParentFile().mkdirs();
        fw = new FileWriter(f);
        velocityEngine.getTemplate(TEMPLATE_TEST_CPP_DIR +"/FunctionH.vm").merge(ctx, fw);
        fw.flush();
        fw.close();
        
	}
}
