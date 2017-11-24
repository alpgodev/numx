/*

TEST VB solution Visual Studio 2003

create:
solution .sln
projects .vbproj

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

public class VB_TEST_VS2003_Generator extends VelocityHelper {

	// paths of velocity models
	public final static String TEMPLATE_TEST_VB_DIR = "com/numx/automatization/generator/test";
	
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
		System.out.println("Generating Visual Studio 2003 solution for VB tests...");
		new VB_TEST_VS2003_Generator().generate(products, outDir);
	}

	public VB_TEST_VS2003_Generator() throws Exception{
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
    	}
	}

	// TEST VB --- create Visual Studio 2003 solution and projects
	public void generateTestByProduct(Product prod, String outTestDir) throws Exception {

		// solution .sln
        File f = new File(outTestDir + "/make/win32/vb-VS2003", prod.get_productShortNameLower() + ".sln");
		f.getParentFile().mkdirs();
        FileWriter fw = new FileWriter(f);
        velocityEngine.getTemplate(TEMPLATE_TEST_VB_DIR +"/Solution-vb-VS2003.vm").merge(ctx, fw);
        fw.flush();
        fw.close(); 

		// projects .vbproj
        f = new File(outTestDir + "/make/win32/vb-VS2003", prod.get_productShortNameLower() + ".vbproj");
		f.getParentFile().mkdirs();
        fw = new FileWriter(f);
        velocityEngine.getTemplate(TEMPLATE_TEST_VB_DIR +"/Project-vb-VS2003.vm").merge(ctx, fw);
        fw.flush();
        fw.close(); 
	}
}
