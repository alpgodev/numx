/*

API C solution Visual Studio 2003 Generator

create:
definition files .def
solution .sln
projects .vcproj

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

public class C_API_VS2003_Generator extends VelocityHelper {

	// paths of velocity models
	public final static String TEMPLATE_MAKE_DIR = "com/numx/automatization/generator/make";
	
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
		System.out.println("Generating API C/C++ solution for Visual Studio 2003 ...");
		new C_API_VS2003_Generator().generate(products, outDir);
	}

	public C_API_VS2003_Generator() throws Exception{
		super();
		VMs.add(TEMPLATE_MAKE_DIR);
        initVelocityEngine();
	}

	public void generate(List products, String outDir) throws Exception {
        ctx = new VelocityContext();
        ctx.put("year", new Integer((new GregorianCalendar()).get(Calendar.YEAR)).toString());

    	for (Iterator itProds = products.iterator(); itProds.hasNext();) {
			Product prod = (Product) itProds.next();
	        ctx.put("prod", prod);
	        
	        String outMakeDir = outDir + "/make/win32";
			
	        generateProjectByProduct(prod, outMakeDir);

			int index = 10;
        	for (Iterator itMods = prod.getModules().iterator(); itMods.hasNext();) {
        		Module mod = (Module) itMods.next();
        		if (mod.is_release()) {
	    	        ctx.put("mod", mod);
					ctx.put("LANG", "c");
					ctx.put("ModProjectGUIDC", "6709D3AC-ED44-43E4-BE68-4DB3D8F4DC"+Integer.toString(index++));
					
					generateProjectByModule(prod, mod, outMakeDir);
	    	        generateMakeByModule(prod, mod, outMakeDir);
        		}
        	}
    	}
	}
	
	// API C --- creation des .def
	public void generateMakeByModule(Product prod, Module mod, String outMakeDir) throws Exception {

		File f = new File(outMakeDir, prod.get_productShortNameLower() + "-" + mod.get_moduleNameLowerScore() + "-c.def");
		f.getParentFile().mkdirs();
        FileWriter fw = new FileWriter(f);
        velocityEngine.getTemplate(TEMPLATE_MAKE_DIR +"/ModuleC.vm").merge(ctx, fw);
        fw.flush();
        fw.close();
        
		String modname = mod.get_moduleNameLowerScore();
        if (modname == "linear-algebra"){
        	modname="linear_algebra";
        }
	}
	
	// API C -- creation de la solution C Visual Studio 2003 (.sln)
	public void generateProjectByProduct(Product prod, String outProjectDir) throws Exception {

		File f = new File(outProjectDir, prod.get_productShortNameLower() + ".sln");
		f.getParentFile().mkdirs();
		FileWriter fw = new FileWriter(f);
		velocityEngine.getTemplate(TEMPLATE_MAKE_DIR +"/Solution-c-VS2003.vm").merge(ctx, fw);
		fw.flush();
		fw.close();	        
       
	}
	
	// API C -- creation des projets C Visual Studio 2003 (.vcproj)
	public void generateProjectByModule(Product prod, Module mod, String outProjectDir) throws Exception {

		File f = new File(outProjectDir, prod.get_productShortNameLower() +"-"+ mod.get_moduleNameLowerScore() +
			"-c.vcproj");
		f.getParentFile().mkdirs();
		FileWriter fw = new FileWriter(f);
		velocityEngine.getTemplate(TEMPLATE_MAKE_DIR +"/Module-c-VS2003.vm").merge(ctx, fw);
		fw.flush();
		fw.close();	  
	}	
}
