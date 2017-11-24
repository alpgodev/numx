/*

 API Java Generator

 create:
  bulidnulx.xml   -> numx/src/api/include/jni
  sources .java   -> numx/src/api/include/jni/com/numx/numx
                  -> numx/src/api/include/jni/com/numx/jee/numx
  sources *jni.c  -> numx/src/api/jni
 

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

public class Java_API_Generator extends VelocityHelper {

	public final static String TEMPLATE_INCLUDE_JNI_DIR = "com/numx/automatization/generator/include";
	public final static String TEMPLATE_SRC_JNI_DIR = "com/numx/automatization/generator/src";
	
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
		System.out.println("Generating API Java ...");
		new Java_API_Generator().generate(products, outDir);
	}

	public Java_API_Generator() throws Exception{
		super();
		System.out.println("	Generating include files .h ...");
		VMs.add(TEMPLATE_INCLUDE_JNI_DIR);
		System.out.println("	Generating source files *jni.c ...");
		VMs.add(TEMPLATE_SRC_JNI_DIR);
        initVelocityEngine();
	}

	public void generate(List products, String outDir) throws Exception {
        ctx = new VelocityContext();
        ctx.put("year", new Integer((new GregorianCalendar()).get(Calendar.YEAR)).toString());

        /* output directory */
        String outIncludeDir = outDir + "/src/api/include/jni";
        
    	for (Iterator itProds = products.iterator(); itProds.hasNext();) {
			Product prod = (Product) itProds.next();
	        ctx.put("prod", prod);
	        
	        generateIncludeByProduct(prod, outIncludeDir);

			int index = 10;
        	for (Iterator itMods = prod.getModules().iterator(); itMods.hasNext();) {
        		Module mod = (Module) itMods.next();
        		if (!mod.is_onlyC() && mod.is_release()) {        			
	    	        ctx.put("mod", mod);
					ctx.put("LANG", "java");
					ctx.put("ModProjectGUIDJAVA", "9F113DAA-866E-4870-B356-EB23C0B57A"+Integer.toString(index++));

					String outSrcDir = outDir + "/src/api/jni/" + mod.get_moduleNameLowerScore();
					
	    	        generateIncludeByModule(prod, mod, outIncludeDir);
	
	            	for (Iterator itFuncs = mod.getFunctions().iterator(); itFuncs.hasNext();) {
	            		Function func = (Function) itFuncs.next();
	            		if (!func.is_onlyC() && func.is_release()) {        		
		        	        ctx.put("func", func);

		        	        generateSrcByFunction(prod, func, outSrcDir);
	            		}
	            	}
        		}
        	}
    	}
	}
	
	/* API Java --- create include files .h */
	public void generateIncludeByProduct(Product prod, String outIncludeDir) throws Exception {

		File f = new File(outIncludeDir + "/com/jee/" + prod.get_productShortNameLower(), prod.get_productShortName() + "ConnectionFactory.java");
		f.getParentFile().mkdirs();
        FileWriter fw = new FileWriter(f);
        velocityEngine.getTemplate(TEMPLATE_INCLUDE_JNI_DIR +"/ConnectionFactoryJEE.vm").merge(ctx, fw);
        fw.flush();
        fw.close();
        
        f = new File(outIncludeDir + "/com/jee/" + prod.get_productShortNameLower(), prod.get_productShortName() + "Connection.java");
		f.getParentFile().mkdirs();
        fw = new FileWriter(f);
        velocityEngine.getTemplate(TEMPLATE_INCLUDE_JNI_DIR +"/ConnectionJEE.vm").merge(ctx, fw);
        fw.flush();
        fw.close();

        f = new File(outIncludeDir + "/com/jee/" + prod.get_productShortNameLower(), prod.get_productShortName() + "ManagedConnectionFactory.java");
		f.getParentFile().mkdirs();
        fw = new FileWriter(f);
        velocityEngine.getTemplate(TEMPLATE_INCLUDE_JNI_DIR +"/ManagedConnectionFactoryJEE.vm").merge(ctx, fw);
        fw.flush();
        fw.close();

        f = new File(outIncludeDir + "/com/jee/" + prod.get_productShortNameLower(), prod.get_productShortName() + "ManagedConnection.java");
		f.getParentFile().mkdirs();
        fw = new FileWriter(f);
        velocityEngine.getTemplate(TEMPLATE_INCLUDE_JNI_DIR +"/ManagedConnectionJEE.vm").merge(ctx, fw);
        fw.flush();
        fw.close();

        f = new File(outIncludeDir + "/com/jee/" + prod.get_productShortNameLower(), prod.get_productShortName() + "ResourceAdapter.java");
		f.getParentFile().mkdirs();
        fw = new FileWriter(f);
        velocityEngine.getTemplate(TEMPLATE_INCLUDE_JNI_DIR +"/ResourceAdapterJEE.vm").merge(ctx, fw);
        fw.flush();
        fw.close();

        f = new File(outIncludeDir + "/jee/META-INF/", prod.get_productShortNameLower() + "_ra.xml");
		f.getParentFile().mkdirs();
        fw = new FileWriter(f);
        velocityEngine.getTemplate(TEMPLATE_INCLUDE_JNI_DIR +"/ra.vm").merge(ctx, fw);
        fw.flush();
        fw.close();

        f = new File(outIncludeDir + "/jee/META-INF/", prod.get_productShortNameLower() + "_jonas-ra.xml");
		f.getParentFile().mkdirs();
        fw = new FileWriter(f);
        velocityEngine.getTemplate(TEMPLATE_INCLUDE_JNI_DIR +"/jonas-ra.vm").merge(ctx, fw);
        fw.flush();
        fw.close();
        
        /* buildnumx.xml */
		f = new File(outIncludeDir, "build" + prod.get_productShortNameLower() + ".xml");
		f.getParentFile().mkdirs();
		fw = new FileWriter(f);
		velocityEngine.getTemplate(TEMPLATE_INCLUDE_JNI_DIR +"/Build.vm").merge(ctx, fw);
		fw.flush();
		fw.close();

	}

    /* API java - */
	public void generateIncludeByModule(Product prod, Module mod, String outIncludeDir) throws Exception {

        /* create NumXModuleName.java -> numx/src/api/include/com/numx/ */
		File f = new File(outIncludeDir + "/com/" + prod.get_productShortNameLower(), prod.get_productShortName() + mod.get_moduleNameNoSpace() + ".java");
		f.getParentFile().mkdirs();
        FileWriter fw = new FileWriter(f);
        velocityEngine.getTemplate(TEMPLATE_INCLUDE_JNI_DIR +"/ModuleJava.vm").merge(ctx, fw);
        fw.flush();
        fw.close();
        
        /* create NumXModuleName.java -> numx/src/api/include/com/jee/numx/ */
        f = new File(outIncludeDir + "/com/jee/" + prod.get_productShortNameLower(), prod.get_productShortName() + mod.get_moduleNameNoSpace() + ".java");
		f.getParentFile().mkdirs();
        fw = new FileWriter(f);
        velocityEngine.getTemplate(TEMPLATE_INCLUDE_JNI_DIR +"/ModuleJEE.vm").merge(ctx, fw);
        fw.flush();
        fw.close();
	}

	/* API Java --- create source files *jni.c -> numx/src/api/jni/ */
	public void generateSrcByFunction(Product prod, Function func, String outSrcDir) throws Exception {

		File f = new File(outSrcDir, func.get_fortranFunction() + "jni.c");
		f.getParentFile().mkdirs();
        FileWriter fw = new FileWriter(f);
        velocityEngine.getTemplate(TEMPLATE_SRC_JNI_DIR +"/FunctionJNI.vm").merge(ctx, fw);
        fw.flush();
        fw.close();
	}
}
