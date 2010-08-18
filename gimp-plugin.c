#include <libgimp/gimp.h>
#include <libgimp/gimpui.h>

typedef struct
{
  gint wsize;
} MyWindowVals;
static gboolean window_dialog (GimpDrawable *drawable,gint32 image_ID);

/* Set up default values for options */
static MyWindowVals bvals =
{
  3  /* Window size default */
};
//Initializing memory
static void
init_mem (guchar ***row,
          gint num_bytes, gint num_rows);
//shuffle rows
static void
shuffle (GimpPixelRgn *rgn_in1,
	GimpPixelRgn *rgn_in2,
	GimpPixelRgn *rgn_in3,
	GimpPixelRgn *rgn_in4,
         guchar      **input1,
	guchar      **input2,
	guchar      **input3,
	guchar      **output,
         gint          x1,
         gint          y1,
         gint          width_m,
         gint          height_m,
         gint          ypos);
//perform calculations
static void
process_window	( guchar **input1,
		  guchar **input2,
		  guchar **input3,
                  guchar **output,
                  gint x1,
		  gint y1,
                  gint width_m,
		  gint height_m,
                  gint channels,
                  gint i);
//query life cycle of plugin
static void query (void);
//runs the plugin
static void run   (const gchar      *name,
                   gint              nparams,
                   const GimpParam  *param,
                   gint             *nreturn_vals,
                   GimpParam       **return_vals);
//reading and writing into image
static void fusion (GimpDrawable     *drawable,gint32 image_ID);

GimpPlugInInfo PLUG_IN_INFO =
{
  NULL,
  NULL,
  query,
  run
};
//global declarations
//Scale of the ms and pan and the window size the user wishes to use
gint scale;
gint wsize;

//for testing
long counter;

MAIN()

static void
query (void)
{
  static GimpParamDef args[] =
  {
    {
      GIMP_PDB_INT32,
      "run-mode",
      "Run mode"
    },
    {
      GIMP_PDB_IMAGE,
      "image",
      "Input image"
    },
    {
      GIMP_PDB_DRAWABLE,
      "drawable",
      "Input drawable"
    }
  };

  gimp_install_procedure (
    "plug-in-fuseimagerow",
    "fuseimagerow",
    "fuseimagerow",
    "DSCE",
    "Copyright PTVV",
    "2008",
    "_fuseimagerow",
    "RGB*, GRAY*",
    GIMP_PLUGIN,
    G_N_ELEMENTS (args), 0,
    args, NULL);

// might need to be changed later the arguments i.e

  gimp_plugin_menu_register ("plug-in-fuseimagerow",
                             "<Image>/Filters/Image Fusion/Fuserow(fast)");
}


static void
run (const gchar      *name,
     gint              nparams,
     const GimpParam  *param,
     gint             *nreturn_vals,
     GimpParam       **return_vals)
{
  static GimpParam  values[1];
  GimpPDBStatusType status = GIMP_PDB_SUCCESS;
  GimpRunMode       run_mode;
  GimpDrawable     *drawable;
  gint32 image_ID;
  /* Setting mandatory output values */
 
  *nreturn_vals = 2;
  *return_vals  = values;

  values[0].type = GIMP_PDB_STATUS;
  values[0].data.d_status = status;

  values[1].type = GIMP_PDB_LAYER;
  values[1].data.d_layer = -1;

  /* Getting run_mode - does not display a dialog if
   * we are in NONINTERACTIVE mode
   */

  run_mode = param[0].data.d_int32;
 image_ID=param[1].data.d_image;
  /*  Get the specified drawable  */
  drawable = gimp_drawable_get (param[2].data.d_drawable);
 switch (run_mode)
    {
    case GIMP_RUN_INTERACTIVE:
      /* Get options last values if needed */
      gimp_get_data ("plug-in-fuseimagerow", &bvals);

      /* Display the dialog */
      if (! window_dialog (drawable,image_ID))
        return;
      break;

    case GIMP_RUN_NONINTERACTIVE:
      if (nparams != 4)
        status = GIMP_PDB_CALLING_ERROR;
      if (status == GIMP_PDB_SUCCESS)
        bvals.wsize = param[3].data.d_int32;
      break;

    case GIMP_RUN_WITH_LAST_VALS:
      /*  Get options last values if needed  */
      gimp_get_data ("plug-in-fuseimagerow", &bvals);
      break;

    default:
      break;
    }
  // progress bar
  gimp_progress_init ("Image Fusion with beta...");


  //image_ID=param[1].data.d_image;

  fusion (drawable,image_ID);


  gimp_displays_flush ();
  gimp_drawable_detach (drawable);
	if (run_mode == GIMP_RUN_INTERACTIVE)
    gimp_set_data ("plug-in-fuseimagerow", &bvals, sizeof (MyWindowVals));
}


static void
fusion (GimpDrawable *drawable,gint32 image_ID)
{
  gint         i,j, channels,m,n=bvals.wsize;
  gint         x1, y1, x2, y2,x11,x12,y11,y12,width_p,height_p,width_m,height_m;
  GimpPixelRgn rgn_in1, rgn_out,rgn_in2,rgn_in3,rgn_in4;
  guchar     **input1,**input2,**input3,**output;
  gint num_layers;
  gint* layers ;
  wsize=n=bvals.wsize;
GimpDrawable* draw0;
GimpDrawable* draw1;
GimpDrawable* draw2;
GimpDrawable* draw3;

	layers=gimp_image_get_layers (image_ID, &num_layers);
	draw3 = gimp_drawable_get (layers[3]);//bottom
	draw2 = gimp_drawable_get (layers[2]);//bottom
        draw1 = gimp_drawable_get (layers[1]);//middle
        draw0 = gimp_drawable_get (layers[0]);//topmost

 gimp_drawable_mask_bounds (draw0->drawable_id,
                            &x1, &y1,
                             &x2, &y2);
gimp_drawable_mask_bounds (draw2->drawable_id,
                            &x11, &y11,
                             &x12, &y12);
width_p  = x12 - x11;
  height_p = y12 - y11;
width_m  = x2 - x1;
  height_m = y2 - y1;
scale=(gint)(width_p/width_m);
m=scale;

wsize=width_m>height_m? height_m:width_m;

  channels = gimp_drawable_bpp (draw0->drawable_id);
  gimp_tile_cache_ntiles (4 * (draw3->width / gimp_tile_width () + 1));

  gimp_pixel_rgn_init (&rgn_in1,
                       draw0,
                       x1, y1,
                       x2 - x1, y2 - y1,
                      FALSE, FALSE);
  gimp_pixel_rgn_init (&rgn_in2,
                       draw1,
                       x1, y1,
                       x2 - x1, y2 - y1,
                      FALSE, FALSE);
  gimp_pixel_rgn_init (&rgn_in3,
                       draw2,
                       x11, y11,
                       x12 - x11, y12 - y11,
                      FALSE, FALSE);

gimp_pixel_rgn_init (&rgn_in4,
                       draw3,
                       x11, y11,
                       x12 - x11, y12 - y11,
                       FALSE, FALSE);
   
gimp_pixel_rgn_init (&rgn_out,
                       draw3,
                       x11, y11,
                       x12 - x11, y12 - y11,
                       TRUE, TRUE);// change this later
  
init_mem (&input1,  width_m * channels, n);
init_mem (&input2,  width_m * channels, n);
init_mem (&input3,  width_p * channels, m*n);
init_mem (&output,  width_p * channels, m*n);

//	g_print("%d\n\n",wsize);

  for (i = 0; i < wsize; i++)
    {
      gimp_pixel_rgn_get_row (&rgn_in1,
                              input1[i],
                              x1, y1 + CLAMP (i, 0, height_m - 1),
                              width_m);
      gimp_pixel_rgn_get_row (&rgn_in2,
                              input2[i],
                              x1, y1 + CLAMP (i, 0, height_m - 1),
                              width_m);
    }
  for (i = 0; i < wsize*scale; i++)
    {
      gimp_pixel_rgn_get_row (&rgn_in3,
                              input3[i],
                              x11, y11 + CLAMP (i, 0, height_p - 1),
                              width_p);

      gimp_pixel_rgn_get_row (&rgn_in4,
                              output[i],
                              x11, y11 + CLAMP (i, 0, height_p - 1),
                              width_p);

    }
  for (i = 0; i < height_m; i+=wsize)
    {
      /* To be done for each tile row */
      process_window (input1,input2,input3,
                   output,
                   x1, y1,
                   width_m, height_m,
                   channels,
                   i);
	for(j=0;j<wsize*scale;j++)
	{
      gimp_pixel_rgn_set_row (&rgn_out,
                              output[j],
                              x11, (scale*i) + j,
                              width_p);
	}
      /* shift tile rows to insert the new one at the end */
      shuffle (&rgn_in1,&rgn_in2,&rgn_in3,&rgn_in4,
               input1,input2,input3,output,
               x1, y1,
               width_m, height_m,
               i+wsize);
      
    }



         gimp_drawable_flush (draw3);
	  gimp_drawable_merge_shadow (draw3->drawable_id, TRUE);
	  gimp_drawable_update (draw3->drawable_id,
                        x11, y11,
                        x12 - x11, y12 - y11);


 for (i = 0; i < wsize; i++)
{
    g_free (input1[i]);
    g_free (input2[i]);
}
for (i = 0; i < wsize*scale; i++)
{
    g_free (input3[i]);
    g_free (output[i]);
}
  g_free (input1);
  g_free (input2);
  g_free (input3);
  g_free (output);


}
	
static void
init_mem (guchar ***row,
          gint num_bytes, gint num_rows)
{
  gint i;

  /* Allocate enough memory for row and outrow */
  *row = g_new (guchar *, num_rows);

  for (i = 0; i < num_rows; i++)
    (*row)[i] = g_new (guchar, num_bytes);
}

static void
process_window	( guchar **input1,
		  guchar **input2,
		  guchar **input3,
                  guchar **output,
                  gint x1,
		  gint y1,
                  gint width_m,
		  gint height_m,
                  gint channels,
                  gint ii)
{
  gint i,j,jj,ssi,ssj,temp;
double sumXY[4],sumX[4],sumY[4],sumX2[4];
double Beta[4];
      gint k,width_p;

width_p=width_m*scale;

wsize=width_m>height_m? height_m:width_m;
//testing
counter=0;
//g_print("%d\n scale = %d\n",wsize,scale);
for (jj = 0; jj < width_m; jj+=wsize)
    {
	

for(k=0;k<4;k++)
	{
        sumXY[k] =0;
	sumX[k] =0;
	sumY[k] =0;
	sumX2[k] =0;
Beta[k]=0;
	}
	for(i=0;i<wsize;i++)
	{
		for(j=0;j<wsize;j++)
		{
			for(k=0;k<4;k++)
			{
	// X= Independent
	// Y= Dependent
        sumXY[k] += (input1[i][channels * CLAMP (j+jj, 0, width_m - 1) + k]*input2[i][channels * CLAMP (j+jj, 0, width_m - 1) + k]);
	sumX[k] += input1[i][channels * CLAMP (j+jj, 0, width_m - 1) + k];
	sumY[k] += input2[i][channels * CLAMP (j+jj, 0, width_m - 1) + k];
	sumX2[k] += (input1[i][channels * CLAMP (j+jj, 0, width_m - 1) + k]*input1[i][channels * CLAMP (j+jj, 0, width_m - 1) + k]);
			}
		}
	}
	for(k=0;k<4;k++)
	{
	//Calcuate Beta for the window 3x3
	if(((wsize*wsize*sumX2[k])-(sumX[k]*sumX[k]))!=0)
       		Beta[k]=((wsize*wsize*sumXY[k])-(sumX[k]*sumY[k]))/((wsize*wsize*sumX2[k])-(sumX[k]*sumX[k]));
	else
	Beta[k]=0;
//	g_print("%f\t",Beta[k]);
	}

	//g_print("%d\t%d\n",x0,y0);
//set the output pixels according to beta
	for(i=0;i<wsize;i++)
	{

		if(width_m>height_m)
       			gimp_progress_update ((gdouble) (i+jj) / (gdouble) width_m);
		else
       			gimp_progress_update ((gdouble) (i+ii) / (gdouble) height_m);

		
	for(j=0;j<wsize;j++)
	{
for(ssi=0;ssi<scale;ssi++){
	for(ssj=0;ssj<scale;ssj++){
		for(k=0;k<4;k++)
		{
// wet
		input3[scale*i+ssi][channels * (CLAMP (scale*(jj+j), 0, width_p - 1)+ssj) + k]= (0.7 * output[scale*i+ssi][channels * (CLAMP (scale*(jj+j), 0, width_p - 1)+ssj) + k]) +( 0.3 * input3[scale*i+ssi][channels * (CLAMP (scale*(jj+j), 0, width_p - 1)+ssj) + k]);
// output	
		temp=input2[i][channels *  CLAMP (j+jj, 0, width_m - 1) + k]+(Beta[k]*(input3[scale*i+ssi][channels * (CLAMP (scale*(jj+j), 0, width_p - 1)+ssj)+k]-input1[i][channels *  CLAMP (j+jj, 0, width_m - 1 + k)]));
if(temp>255)
				output[scale*i+ssi][channels * (CLAMP (scale*(jj+j), 0, width_p - 1)+ssj) + k]=255;
else if(temp<1)
				output[scale*i+ssi][channels * (CLAMP (scale*(jj+j), 0, width_p - 1)+ssj) + k]=0;
else
				output[scale*i+ssi][channels * (CLAMP (scale*(jj+j), 0, width_p - 1)+ssj) + k]=temp;
	
	


}


//pixel discoloration testing
		if((output[scale*i+ssi][channels * (CLAMP (scale*(jj+j), 0, width_p - 1)+ssj) + 0]-input2[i][channels *  CLAMP (j+jj, 0, width_m - 1) + 0])>220 ||
		(output[scale*i+ssi][channels * (CLAMP (scale*(jj+j), 0, width_p - 1)+ssj) + 1]-input2[i][channels *  CLAMP (j+jj, 0, width_m - 1) + 1])>220 ||
	(output[scale*i+ssi][channels * (CLAMP (scale*(jj+j), 0, width_p - 1)+ssj) + 2]-input2[i][channels *  CLAMP (j+jj, 0, width_m - 1) + 2])>220)
		{
			if((output[scale*i+ssi][channels * (CLAMP (scale*(jj+j), 0, width_p - 1)+ssj) + 0]-input2[i][channels *  CLAMP (j+jj, 0, width_m - 1) + 0])>220)
			{
		
output[scale*i+ssi][channels * (CLAMP (scale*(jj+j), 0, width_p - 1)+ssj) + 0] = (0.7 * input3[scale*i+ssi][channels * (CLAMP (scale*(jj+j), 0, width_p - 1)+ssj) + 0]) + (0.3 * input2[i][channels *  CLAMP (j+jj, 0, width_m - 1) + 0]);
			}

			if((output[scale*i+ssi][channels * (CLAMP (scale*(jj+j), 0, width_p - 1)+ssj) + 1]-input2[i][channels *  CLAMP (j+jj, 0, width_m - 1) + 1])>220)
			{
	output[scale*i+ssi][channels * (CLAMP (scale*(jj+j), 0, width_p - 1)+ssj) + 1] = (0.7 * input3[scale*i+ssi][channels * (CLAMP (scale*(jj+j), 0, width_p - 1)+ssj) + 1]) + (0.3 * input2[i][channels *  CLAMP (j+jj, 0, width_m - 1) + 1]);

			}

			if((output[scale*i+ssi][channels * (CLAMP (scale*(jj+j), 0, width_p - 1)+ssj) + 2]-input2[i][channels *  CLAMP (j+jj, 0, width_m - 1) + 2])>220)
			{
	output[scale*i+ssi][channels * (CLAMP (scale*(jj+j), 0, width_p - 1)+ssj) + 2] = (0.7 * input3[scale*i+ssi][channels * (CLAMP (scale*(jj+j), 0, width_p - 1)+ssj) + 2]) + (0.3 * input2[i][channels *  CLAMP (j+jj, 0, width_m - 1) + 2]);

			}
		}

				}
			}	

	}
	}

    }
}

static void
shuffle (GimpPixelRgn *rgn_in1,
	GimpPixelRgn *rgn_in2,
	GimpPixelRgn *rgn_in3,
	GimpPixelRgn *rgn_in4,
         guchar      **input1,
	guchar      **input2,
	guchar      **input3,
	guchar      **output,
         gint          x1,
         gint          y1,
         gint          width_m,
         gint          height_m,
         gint          ypos)
{
  gint    i;
  gint    width_p=width_m*scale, height_p=height_m*scale;

  for (i = 0; i < wsize; i++)
    {
      gimp_pixel_rgn_get_row (rgn_in1,
                              input1[i],
                              x1, MIN (ypos + i + y1, y1 + height_m - 1),
                              width_m);
      gimp_pixel_rgn_get_row (rgn_in2,
                              input2[i],
                              x1, MIN (ypos + i + y1, y1 + height_m - 1),
                              width_m);
    }
  
  for (i = 0; i < wsize*scale; i++)
    {
      gimp_pixel_rgn_get_row (rgn_in3,
                              input3[i],
                              x1, MIN ((scale*ypos) + i + y1, y1 + height_p - 1),
                              width_p);
      gimp_pixel_rgn_get_row (rgn_in4,
                              output[i],
                              x1, MIN ((scale*ypos) + i + y1, y1 + height_p - 1),
                              width_p);

    }

}
static gboolean
window_dialog (GimpDrawable *drawable,gint32 image_ID)
{
  GtkWidget *dialog;
  GtkWidget *main_vbox;
  GtkWidget *main_hbox;
  GtkWidget *frame;
  GtkWidget *radius_label;
  GtkWidget *alignment;
  GtkWidget *spinbutton;
  GtkObject *spinbutton_adj;
  GtkWidget *frame_label;
  GimpDrawable* draw1;
  gint num_layers;
  gint* layers ;
  gint x1,x2,y1,y2;
  gboolean   run;
	layers=gimp_image_get_layers (image_ID, &num_layers);
        draw1 = gimp_drawable_get (layers[1]);//middle
	 gimp_drawable_mask_bounds (draw1->drawable_id,
                            &x1, &y1,
                             &x2, &y2);

  gimp_ui_init ("fuseimagerow", FALSE);

  dialog = gimp_dialog_new ("Fuseimagerow", "fuseimagerow",
                            NULL, 0,
                            gimp_standard_help_func, "plug-in-fuseimagerow",

                            GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                            GTK_STOCK_OK,     GTK_RESPONSE_OK,

                            NULL);

  main_vbox = gtk_vbox_new (FALSE, 6);
  gtk_container_add (GTK_CONTAINER (GTK_DIALOG (dialog)->vbox), main_vbox);
  gtk_widget_show (main_vbox);

  frame = gtk_frame_new (NULL);
  gtk_widget_show (frame);
  gtk_box_pack_start (GTK_BOX (main_vbox), frame, TRUE, TRUE, 0);
  gtk_container_set_border_width (GTK_CONTAINER (frame), 6);

  alignment = gtk_alignment_new (0.5, 0.5, 1, 1);
  gtk_widget_show (alignment);
  gtk_container_add (GTK_CONTAINER (frame), alignment);
  gtk_alignment_set_padding (GTK_ALIGNMENT (alignment), 6, 6, 6, 6);

  main_hbox = gtk_hbox_new (FALSE, 0);
  gtk_widget_show (main_hbox);
  gtk_container_add (GTK_CONTAINER (alignment), main_hbox);

  radius_label = gtk_label_new_with_mnemonic ("_Wsize:");
  gtk_widget_show (radius_label);
  gtk_box_pack_start (GTK_BOX (main_hbox), radius_label, FALSE, FALSE, 6);
  gtk_label_set_justify (GTK_LABEL (radius_label), GTK_JUSTIFY_RIGHT);

  spinbutton_adj = gtk_adjustment_new (3, 1, (x2>y2?y2:x2), 1, 5, 5);


  spinbutton = gtk_spin_button_new (GTK_ADJUSTMENT (spinbutton_adj), 1, 0);
  gtk_widget_show (spinbutton);
  gtk_box_pack_start (GTK_BOX (main_hbox), spinbutton, FALSE, FALSE, 6);
  gtk_spin_button_set_numeric (GTK_SPIN_BUTTON (spinbutton), TRUE);

  frame_label = gtk_label_new ("<b>Modify Window Size</b>");
  gtk_widget_show (frame_label);
  gtk_frame_set_label_widget (GTK_FRAME (frame), frame_label);
  gtk_label_set_use_markup (GTK_LABEL (frame_label), TRUE);

  g_signal_connect (spinbutton_adj, "value_changed",
                    G_CALLBACK (gimp_int_adjustment_update),
                    &bvals.wsize);
  gtk_widget_show (dialog);

  run = (gimp_dialog_run (GIMP_DIALOG (dialog)) == GTK_RESPONSE_OK);

  gtk_widget_destroy (dialog);
bvals.wsize=x2>y2?y2:x2;
  return run;
}


