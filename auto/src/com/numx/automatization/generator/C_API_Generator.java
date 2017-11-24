/***********************
 * API C Generator     *
 *                     *
 * create:             *
 *  includes files .h  *
 *  sources .c         *
 ***********************/

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

public class C_API_Generator extends VelocityHelper {

	// paths of velocity models
	public final static String TEMPLATE_INCLUDE_C_DIR = "com/numx/automatization/generator/include";
	public final static String TEMPLATE_SRC_C_DIR = "com/numx/automatization/generator/src";
	
	private Context ctx;
	
	/**
	 * @param args
	 * args[0] is the path of descriptor files
	 * args[1] is the path of root of numX
	 */
	public static void main(String[] args) throws Exception{
		//reading descripor files and getting a list of elements of type Product
		XMLParser parser = new XMLParser();
		List products = parser.getProducts(args[0]);

		//generating files
		String outDir = args[1];
		System.out.println("Generating API C/C++ ...");
		new C_API_Generator().generate(products, outDir);
	}

	public C_API_Generator() throws Exception{
		super();
		System.out.println("	Generating include files .h ...");
		VMs.add(TEMPLATE_INCLUDE_C_DIR);
		System.out.println("	Generating source files .c ...");
		VMs.add(TEMPLATE_SRC_C_DIR);
        initVelocityEngine();
	}

	public void generate(List products, String outDir) throws Exception {
        ctx = new VelocityContext();
        ctx.put("year", new Integer((new GregorianCalendar()).get(Calendar.YEAR)).toString());

        String outIncludeDir = outDir + "/lib/include";
    	for (Iterator itProds = products.iterator(); itProds.hasNext();) {
			Product prod = (Product) itProds.next();
	        ctx.put("prod", prod);

			int index = 10;
        	for (Iterator itMods = prod.getModules().iterator(); itMods.hasNext();) {
        		Module mod = (Module) itMods.next();
        		if (mod.is_release()) {
	    	        ctx.put("mod", mod);
					ctx.put("LANG", "c");
					ctx.put("ModProjectGUIDC", "6709D3AC-ED44-43E4-BE68-4DB3D8F4DC"+Integer.toString(index++));

					String outSrcDir = outDir + "/src/api/c/" + mod.get_moduleNameLowerScore();
					
	    	        generateIncludeByModule(prod, mod, outIncludeDir);
	
	            	for (Iterator itFuncs = mod.getFunctions().iterator(); itFuncs.hasNext();) {
	            		Function func = (Function) itFuncs.next();
	            		if (func.is_release()) {        		
		        	        ctx.put("func", func);
		        	        generateSrcByFunction(prod, func, outSrcDir);
	            		}
	            	}
        		}
        	}
    	}
	}
	
	/* API C --- create include files .h */
	public void generateIncludeByModule(Product prod, Module mod, String outIncludeDir) throws Exception {

		File f = new File(outIncludeDir, prod.get_productShortName() + mod.get_moduleNameNoSpace() + ".h");
		f.getParentFile().mkdirs();
        FileWriter fw = new FileWriter(f);
        velocityEngine.getTemplate(TEMPLATE_INCLUDE_C_DIR +"/ModuleC.vm").merge(ctx, fw);
        fw.flush();
        fw.close();

	}

	/* API C --- create source files .c */
	public void generateSrcByFunction(Product prod, Function func, String outSrcDir) throws Exception {

		File f = new File(outSrcDir, func.get_fortranFunction() + "c.c");
		f.getParentFile().mkdirs();
        FileWriter fw = new FileWriter(f);
        velocityEngine.getTemplate(TEMPLATE_SRC_C_DIR +"/FunctionC.vm").merge(ctx, fw);
        fw.flush();
        fw.close();
	}
}
