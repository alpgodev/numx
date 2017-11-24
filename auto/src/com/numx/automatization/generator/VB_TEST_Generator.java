/*

TEST VB Generator

create:
sources .vb

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

public class VB_TEST_Generator extends VelocityHelper {

	public final static String TEMPLATE_TEST_VB_DIR = "com/numx/automatization/generator/test";
	
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
		System.out.println("Generating tests for Visual Basic ...");
		new VB_TEST_Generator().generate(products, outDir);
	}

	public VB_TEST_Generator() throws Exception{
		super();
		VMs.add(TEMPLATE_TEST_VB_DIR);
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
					ctx.put("LANG", "vb");
					ctx.put("ModProjectGUIDVB", "9C1E0C96-F3BA-4385-A62F-4F85D6A86B"+Integer.toString(index++));

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

	// TESTS VB --- creation des fichiers product pour les tests
	public void generateTestByProduct(Product prod, String outTestDir) throws Exception {

        File f = new File(outTestDir + "/src/vb-net", "Test" + prod.get_productShortName() + ".vb");
		f.getParentFile().mkdirs();
        FileWriter fw = new FileWriter(f);
        velocityEngine.getTemplate(TEMPLATE_TEST_VB_DIR +"/Product-vb-VS2003.vm").merge(ctx, fw);
        fw.flush();
        fw.close();

	}

	// TESTS VB --- creation des fichiers modules .vb pour les tests
	public void generateTestByModule(Product prod, Module mod, String outTestDir) throws Exception {
		
		File f = new File(outTestDir + "/src/vb-net", "Test" + prod.get_productShortName() + mod.get_moduleNameNoSpace() + ".vb");
		f.getParentFile().mkdirs();
        FileWriter fw = new FileWriter(f);
        velocityEngine.getTemplate(TEMPLATE_TEST_VB_DIR +"/Module-vb-VS2003.vm").merge(ctx, fw);
        fw.flush();
        fw.close();

	}

	// TESTS VB -- creation des fichiers fonctions .vb pour les tests
	public void generateTestByFunction(Function func, String outTestDir) throws Exception {
	
        File f = new File(outTestDir + "/src/vb-net/", "Test" + func.get_functionName() + ".vb");
		f.getParentFile().mkdirs();
        FileWriter fw = new FileWriter(f);
        velocityEngine.getTemplate(TEMPLATE_TEST_VB_DIR +"/Function-vb-VS2003.vm").merge(ctx, fw);
        fw.flush();
        fw.close();
        
	}
}
