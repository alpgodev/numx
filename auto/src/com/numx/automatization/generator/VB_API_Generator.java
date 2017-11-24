/*

 API VB Generator

 create:
  documentation HTLM
  includes files .bas .vb
  sources .c

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

public class VB_API_Generator extends VelocityHelper {

	public final static String TEMPLATE_HTML_DOC_DIR = "com/numx/automatization/generator/doc";
	public final static String TEMPLATE_INCLUDE_VB_DIR = "com/numx/automatization/generator/include";
	public final static String TEMPLATE_SRC_C_DIR = "com/numx/automatization/generator/src";
	
	private Context ctx;
	
	/**
	 * @param args
	 * args[0] is the path of descriptor files
	 * args[1] is the path of root of numX
	 */
	public static void main(String[] args) throws Exception 
	{
		//reading descripor files and getting a list of elements of type Product
		XMLParser parser = new XMLParser();
		List products = parser.getProducts(args[0]);
		
		//generating files
		String outDir = args[1];
		System.out.println("Generating API VB ...");
		new VB_API_Generator().generate(products, outDir);
	}

	public VB_API_Generator() throws Exception{
		super();
		System.out.println("	Generating documentation .htlm ...");
		VMs.add(TEMPLATE_HTML_DOC_DIR);
		System.out.println("	Generating include files .h ...");
		VMs.add(TEMPLATE_INCLUDE_VB_DIR);
		System.out.println("	Generating sources files .c ...");
		VMs.add(TEMPLATE_SRC_C_DIR);
        initVelocityEngine();
	}

	public void generate(List products, String outDir) throws Exception {
        ctx = new VelocityContext();
        ctx.put("year", new Integer((new GregorianCalendar()).get(Calendar.YEAR)).toString());

        String outIncludeDir = outDir + "/release/api";
        
    	for (Iterator itProds = products.iterator(); itProds.hasNext();) {
			Product prod = (Product) itProds.next();
	        ctx.put("prod", prod);
	        
	        String outHtmlDocDir = outDir + "/doc/manual/html";
	        String outMakeDir = outDir + "/make/win32";

	        generateHtmlDocByProduct(prod, outHtmlDocDir);
	        generateIncludeByProduct(prod, outIncludeDir);

			int index = 10;
        	for (Iterator itMods = prod.getModules().iterator(); itMods.hasNext();) {
        		Module mod = (Module) itMods.next();
        		if (!mod.is_onlyC() && mod.is_release()) {        			
	    	        ctx.put("mod", mod);
					ctx.put("LANG", "vb");
					ctx.put("ModProjectGUIDVB", "9C1E0C96-F3BA-4385-A62F-4F85D6A86B"+Integer.toString(index++));

					String outSrcDir = outDir + "/src/api/c/" + mod.get_moduleNameLowerScore();

	    	        generateHtmlDocByModule(mod, outHtmlDocDir);
	    	        generateIncludeByModule(prod, mod, outIncludeDir);
	
	            	for (Iterator itFuncs = mod.getFunctions().iterator(); itFuncs.hasNext();) {
	            		Function func = (Function) itFuncs.next();
	            		if (!func.is_onlyC() && func.is_release()) {        		
		        	        ctx.put("func", func);

		        	        generateHtmlDocByFunction(func, outHtmlDocDir);
		        	        generateSrcByFunction(prod, func, outSrcDir);
	            		}
	            	}
        		}
        	}
    	}
	}

	// create documentation .hhp .hhc
	public void generateHtmlDocByProduct(Product prod, String outHtmlDocDir) throws Exception {

        File f = new File(outHtmlDocDir, "FunctionsDescription.html");
		f.getParentFile().mkdirs();
        FileWriter fw = new FileWriter(f);
        velocityEngine.getTemplate(TEMPLATE_HTML_DOC_DIR +"/Product.vm").merge(ctx, fw);
        fw.flush();
        fw.close();	        

        f = new File(outHtmlDocDir, prod.get_productShortNameLower() + ".hhp");
		f.getParentFile().mkdirs();
        fw = new FileWriter(f);
        velocityEngine.getTemplate(TEMPLATE_HTML_DOC_DIR +"/Project.vm").merge(ctx, fw);
        fw.flush();
        fw.close();	        

        f = new File(outHtmlDocDir, prod.get_productShortNameLower() + ".hhc");
		f.getParentFile().mkdirs();
        fw = new FileWriter(f);
        velocityEngine.getTemplate(TEMPLATE_HTML_DOC_DIR +"/Content.vm").merge(ctx, fw);
        fw.flush();
        fw.close();	        
	}

	// create documentation .html
	public void generateHtmlDocByModule(Module mod, String outHtmlDocDir) throws Exception {

		File f = new File(outHtmlDocDir, mod.get_moduleNameNoSpace() + ".html");
		f.getParentFile().mkdirs();
        FileWriter fw = new FileWriter(f);
        velocityEngine.getTemplate(TEMPLATE_HTML_DOC_DIR +"/Module.vm").merge(ctx, fw);
        fw.flush();
        fw.close();
	}
	public void generateHtmlDocByFunction(Function func, String outHtmlDocDir) throws Exception {

		File f = new File(outHtmlDocDir, func.get_functionName() + ".html");
		f.getParentFile().mkdirs();
        FileWriter fw = new FileWriter(f);
        velocityEngine.getTemplate(TEMPLATE_HTML_DOC_DIR +"/Function.vm").merge(ctx, fw);
        fw.flush();
        fw.close();
	}

	// Error Codes - create VB description files (.bas)  
	public void generateIncludeByProduct(Product prod, String outIncludeDir) throws Exception {

		File f = new File(outIncludeDir + "/vb-6", "NumxErrorCodes.bas");
		f.getParentFile().mkdirs();
        FileWriter fw = new FileWriter(f);
        velocityEngine.getTemplate(TEMPLATE_INCLUDE_VB_DIR +"/ErrorCodesVB.vm").merge(ctx, fw);
        fw.flush();
        fw.close();

	}

	// Product Module - create VB description files .bas .vb
	public void generateIncludeByModule(Product prod, Module mod, String outIncludeDir) throws Exception {

		File f = new File(outIncludeDir + "/vb-6", prod.get_productShortName() + mod.get_moduleNameNoSpace() + ".bas");
		f.getParentFile().mkdirs();
        FileWriter fw = new FileWriter(f);
        velocityEngine.getTemplate(TEMPLATE_INCLUDE_VB_DIR +"/ModuleVB.vm").merge(ctx, fw);
        fw.flush();
        fw.close();

		f = new File(outIncludeDir + "/vb-net", prod.get_productShortName() + mod.get_moduleNameNoSpace() + ".vb");
		f.getParentFile().mkdirs();
        fw = new FileWriter(f);
        velocityEngine.getTemplate(TEMPLATE_INCLUDE_VB_DIR +"/ModuleVBNET.vm").merge(ctx, fw);
        fw.flush();
        fw.close();

	}

	// API VB --- create source files .c
	public void generateSrcByFunction(Product prod, Function func, String outSrcDir) throws Exception {

		File f = new File(outSrcDir, func.get_fortranFunction() + "c.c");
		f.getParentFile().mkdirs();
        FileWriter fw = new FileWriter(f);
        velocityEngine.getTemplate(TEMPLATE_SRC_C_DIR +"/FunctionC.vm").merge(ctx, fw);
        fw.flush();
        fw.close();
	}
}
