/*

TEST C solution Visual Studio 2003

create:
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

public class C_TEST_VS2003_Generator extends VelocityHelper {

	// paths of velocity models
	public final static String TEMPLATE_TEST_CPP_DIR = "com/numx/automatization/generator/test";
	
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
		System.out.println("Generating Visual Studio 2003 solution for C tests...");
		new C_TEST_VS2003_Generator().generate(products, outDir);
	}

	public C_TEST_VS2003_Generator() throws Exception{
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
    	}
	}
	
	// TEST C --- create Visual Studio 2003 solution and projects
	public void generateTestByProduct(Product prod, String outTestDir) throws Exception {

		// solution C Visual Studio 2003 (.sln)
        File f = new File(outTestDir + "/make/win32/c-VS2003", prod.get_productShortNameLower() + ".sln");
		f.getParentFile().mkdirs();
        FileWriter fw = new FileWriter(f);
        velocityEngine.getTemplate(TEMPLATE_TEST_CPP_DIR +"/Solution-c-VS2003.vm").merge(ctx, fw);
        fw.flush();
        fw.close(); 

		// projets C Visual Studio 2003 (.vcproj)
        f = new File(outTestDir + "/make/win32/c-VS2003", prod.get_productShortNameLower() + ".vcproj");
		f.getParentFile().mkdirs();
        fw = new FileWriter(f);
        velocityEngine.getTemplate(TEMPLATE_TEST_CPP_DIR +"/Project-c-VS2003.vm").merge(ctx, fw);
        fw.flush();
        fw.close();

	}
}
