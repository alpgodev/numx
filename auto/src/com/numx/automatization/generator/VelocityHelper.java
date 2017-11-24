/**
 * Copyright (C) 2012
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */
package com.numx.automatization.generator;

import org.apache.velocity.app.Velocity;
import org.apache.velocity.app.VelocityEngine;

import java.util.ArrayList;
import java.util.List;

public class VelocityHelper {

    /**
     * Velocity constant.
     */
    private static final String CLASS_RESOURCE_LOADER_DESCRIPTION_PROP = "class.resource.loader.description";

    /**
     * Velocity constant.
     */
    private static final String CLASS_RESOURCE_LOADER_DESCRIPTION_VALUE = "Velocity Classpath Resource Loader";

    /**
     * Velocity constant.
     */
    private static final String CLASS_RESOURCE_LOADER_CLASS_PROP = "class.resource.loader.class";

    /**
     * Velocity constant.
     */
    private static final String CLASS_RESOURCE_LOADER_CLASS_VALUE = "org.apache.velocity.runtime.resource.loader.ClasspathResourceLoader";

    public VelocityEngine velocityEngine = null;
    
    public List VMs = new ArrayList();
    
    public void initVelocityEngine() throws Exception {
        if (velocityEngine == null) {
            velocityEngine = new VelocityEngine();
            velocityEngine.setProperty(Velocity.RESOURCE_LOADER, "class");
            velocityEngine.setProperty(CLASS_RESOURCE_LOADER_DESCRIPTION_PROP,
                    CLASS_RESOURCE_LOADER_DESCRIPTION_VALUE);
            velocityEngine.setProperty(CLASS_RESOURCE_LOADER_CLASS_PROP,
                    CLASS_RESOURCE_LOADER_CLASS_VALUE);

            if (VMs.size() > 0) {
                StringBuffer strlibs = new StringBuffer();
                String sep = "";
                for (int i = 0; i < VMs.size(); i++) {
                    strlibs.append(sep);
                    sep = ",";
                    strlibs.append(VMs.get(i));
                }
                velocityEngine.setProperty(Velocity.VM_LIBRARY, strlibs.toString());
            }
            velocityEngine.init();
        }
    }
}
