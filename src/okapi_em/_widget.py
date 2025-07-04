'''
Copyright 2022 Rosalind Franklin Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
'''

'''
Main napari widget with data selector and tabs with the tools
'''

#from qtpy import QtWidgets
from qtpy.QtWidgets import QWidget, QTabWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QGroupBox, QGridLayout
from qtpy.QtWidgets import QRadioButton, QCheckBox,QPushButton, QLabel,QSpinBox, QTableWidget, QTableWidgetItem
from qtpy.QtWidgets import QDoubleSpinBox, QComboBox
from qtpy.QtCore import Qt

import napari
from napari.types import ImageData, LayerDataTuple
from napari.utils import progress

#from magicgui.widgets import Table #Note, this is not documented in magicgui

from . import filters
from . import measure_charging

import numpy as np

# from magicgui import magicgui
# from magicgui import magic_factory

import chafer
from . import quoll_wrapper as quoll

from . import slice_alignment

class MainQWidget(QWidget):
    # Sets the main window of the widget
    # laying out all controls and doing relevant connections between elements and code

    curData = None
    curDataSlice = None
    myFFTFilterDirSliceRes = None


    def __init__(self, napari_viewer: 'napari.viewer.Viewer'):
        super().__init__()
        self.viewer = napari_viewer #Reference to the viewer that will be needed thoughout the class

        #print(f"napari_viewer: {napari_viewer}")

        self.setLayout(QVBoxLayout())
        
        # Set default data source button parent
        self.w_setselect = widget_ImageSetCurrentSelect(napari_viewer, "Data source")

        #Creates the user interface 
        self.layout().addWidget(self.w_setselect)

        #Put a tab widget that will help select operations
        self.tabwidget = QTabWidget()
        self.layout().addWidget(self.tabwidget)
        # and set its layout as vertical
        #self.t0_vbox0.addWidget(self.tabwidget)

        #Create all tabs here
        self.ui_do_tab_SliceAlignment()
        self.ui_do_tab_ChargeSuppression()
        self.ui_do_tab_ResMeas()
        self.ui_do_tab_Filters()
    
    def ui_do_tab_ChargeSuppression(self):
        # /tabChargeSuppr \ t0_vbox0
        self.tabChargeSuppr = QWidget()
        self.tabwidget.addTab(self.tabChargeSuppr, "C-Suppression")
        # Tab with charge suppression filters
        # Contains FFT directional filter plus others

        t0_vbox0 = QVBoxLayout() #Items organised vertically
        self.tabChargeSuppr.setLayout(t0_vbox0)

        # Filter: Directional band-pass
        # Radiobutton followed by a box with controls for this filter
        self.rbDirFFT = QRadioButton("Directional band-pass FFT filter")
        self.rbDirFFT.setChecked(False)
        t0_vbox0.addWidget(self.rbDirFFT)
        #Options for this filter are contained in the next box grpBox1
        grpBox1 = QGroupBox("") #Box, no title, radiobutton does it
        #self.tabChargeSuppr.addWidget(grpBox1)
        t0_vbox0.addWidget(grpBox1)

        vbox2 = QVBoxLayout()
        grpBox1.setLayout(vbox2)
        label1 = QLabel("Select cutoff frequencies")
        vbox2.addWidget(label1)
        t0_hbox = QHBoxLayout() #Horizontal layout to set low and high freq of filters
        vbox2.addLayout( t0_hbox)
        label_low = QLabel("low:")
        t0_hbox.addWidget(label_low)
        self.freqselector_low = QSpinBox()
        self.freqselector_low.setMinimum(0)
        self.freqselector_low.setMaximum(100)
        self.freqselector_low.setValue(0)
        self.freqselector_low.setSingleStep(1)
        t0_hbox.addWidget(self.freqselector_low)
        self.freqselector_low.valueChanged.connect(self.freqselector_low_valueChanged)
        self.label_high = QLabel("high:")
        t0_hbox.addWidget(self.label_high)
        self.freqselector_high = QSpinBox()
        self.freqselector_high.setMinimum(0)
        self.freqselector_high.setMaximum(100)
        self.freqselector_high.setValue(5)
        self.freqselector_high.setSingleStep(1)
        t0_hbox.addWidget(self.freqselector_high)
        self.freqselector_high.valueChanged.connect(self.freqselector_high_valueChanged)

        self.rbChafer = QRadioButton("Chafer")
        self.rbChafer.setChecked(True)
        t0_vbox0.addWidget(self.rbChafer)
        self.grpBoxChafer = QGroupBox("")
        t0_vbox0.addWidget(self.grpBoxChafer) #Add Chafer UI to this grpBoxChafer
        self.grpBoxChafer_layout = QVBoxLayout()
        self.grpBoxChafer.setLayout(self.grpBoxChafer_layout)

        #label data selector
        self.w_setselectlabel= widget_ImageSetCurrentSelect(self.viewer, "Labels data source")
        self.grpBoxChafer_layout.addWidget(self.w_setselectlabel)

        #need parameters
        #nlinesaverage=20,
        #data_fit_max_length_edge_px = 700,
        #data_fit_min_length_px = 50

        qgrid2_layout = QGridLayout()
        self.grpBoxChafer_layout.addLayout(qgrid2_layout)

        qgrid2_layout.addWidget( QLabel("nlinesaverage") ,0,0)
        self.chafer_widg_nlines = QSpinBox()
        self.chafer_widg_nlines.setMinimum(1)
        self.chafer_widg_nlines.setMaximum(65536)
        self.chafer_widg_nlines.setValue(20)
        qgrid2_layout.addWidget( self.chafer_widg_nlines ,0,1)

        qgrid2_layout.addWidget( QLabel("fit length min px") ,1,0)
        self.chafer_widg_minlength = QSpinBox()
        self.chafer_widg_minlength.setMinimum(1)
        self.chafer_widg_minlength.setMaximum(65536)
        self.chafer_widg_minlength.setValue(50)
        qgrid2_layout.addWidget( self.chafer_widg_minlength ,1,1)

        qgrid2_layout.addWidget( QLabel("fit length max px") ,2,0)
        self.chafer_widg_maxlength = QSpinBox()
        self.chafer_widg_maxlength.setMinimum(1)
        self.chafer_widg_maxlength.setMaximum(65536)
        self.chafer_widg_maxlength.setValue(700)
        qgrid2_layout.addWidget( self.chafer_widg_maxlength ,2,1)

        #Need select source file,
        #Need user to select labels file

        self.btnFFTDirFilterSlice = QPushButton("Filter (slice)")
        t0_vbox0.addWidget(self.btnFFTDirFilterSlice)
        self.btnFFTDirFilterSlice.clicked.connect(self.doFFTDirFilterSlice)

        self.btnApply = QPushButton("Filter")
        t0_vbox0.addWidget(self.btnApply)
        self.btnApply.clicked.connect(self.doFFTDirFilterWhole) #Signal

    def ui_do_tab_SliceAlignment(self):
        # /tabSliceAlignemnt \
        self.tabSliceAlignment = QWidget()
        self.tabwidget.addTab(self.tabSliceAlignment, "Slice Alignment")
        tabSliceAlignment_l= QVBoxLayout() #Items organised vertically
        self.tabSliceAlignment.setLayout(tabSliceAlignment_l)
        
        self.chkbxTranslate= QCheckBox("Translation")
        self.chkbxTranslate.setChecked(True)
        tabSliceAlignment_l.addWidget(self.chkbxTranslate)
        self.chkbxLinear= QCheckBox("Linear all free (ov. below)")
        tabSliceAlignment_l.addWidget(self.chkbxLinear)

        hl_Rot = QHBoxLayout()
        tabSliceAlignment_l.addLayout(hl_Rot)
        self.chkbxRotation = QCheckBox("rotation (ov. shear)")
        hl_Rot.addWidget(self.chkbxRotation)
        self.chkbxShearX = QCheckBox("shear X")
        self.chkbxShearX.setChecked(True)
        hl_Rot.addWidget(self.chkbxShearX)
        self.chkbxShearY = QCheckBox("shear Y")
        hl_Rot.addWidget(self.chkbxShearY)

        hl_Scaling = QHBoxLayout()
        tabSliceAlignment_l.addLayout(hl_Scaling)
        self.chkbxScaling = QCheckBox("scaling (ov. stretch)")
        hl_Scaling.addWidget(self.chkbxScaling)
        self.chkbxStretchX = QCheckBox("stretch x")
        hl_Scaling.addWidget(self.chkbxStretchX)
        self.chkbxStretchY = QCheckBox("stretch y")
        hl_Scaling.addWidget(self.chkbxStretchY)
        self.chkbxStretchY.setChecked(True)

        #parameters to control slice alignment feature matching
        qgrid_adj_params = QGridLayout()
        tabSliceAlignment_l.addLayout(qgrid_adj_params)
        label_slice_align0 = QLabel("Ratio test:")
        qgrid_adj_params.addWidget(label_slice_align0,0,0)
        self.spb_SliceAlign_RatioTest=QDoubleSpinBox()
        self.spb_SliceAlign_RatioTest.setMinimum(0.1)
        self.spb_SliceAlign_RatioTest.setMaximum(1.0)
        self.spb_SliceAlign_RatioTest.setDecimals(2)
        self.spb_SliceAlign_RatioTest.setSingleStep(0.1)
        self.spb_SliceAlign_RatioTest.setValue(0.80)
        self.spb_SliceAlign_RatioTest.setToolTip("Ratio test filters matching features based in the ratio to second match score.\nHigher value gets more matches but lower match confidence.\nValues from 0.1 to 1.0.")
        qgrid_adj_params.addWidget(self.spb_SliceAlign_RatioTest, 0,1)

        label_d2dmedian_ratio = QLabel("D2Dmedian ratio:")
        qgrid_adj_params.addWidget(label_d2dmedian_ratio,1,0)
        self.spb_SliceAlign_D2Dmedian_Ratio=QDoubleSpinBox()
        self.spb_SliceAlign_D2Dmedian_Ratio.setMinimum(0.01)
        self.spb_SliceAlign_D2Dmedian_Ratio.setMaximum(10.0)
        self.spb_SliceAlign_D2Dmedian_Ratio.setDecimals(2)
        self.spb_SliceAlign_D2Dmedian_Ratio.setSingleStep(0.01)
        self.spb_SliceAlign_D2Dmedian_Ratio.setValue(0.05)
        self.spb_SliceAlign_D2Dmedian_Ratio.setToolTip("Distance to median ratio filter.\nHigher value gets more matches but less precision.\nValues from 0.001 to 10.0. default 0.05")
        qgrid_adj_params.addWidget(self.spb_SliceAlign_D2Dmedian_Ratio, 1,1)

        self.btn_testSA=QPushButton("Test/preview slice alignment")
        qgrid_adj_params.addWidget(self.btn_testSA, 2,0)
        self.btn_testSA.clicked.connect(self.btn_testSA_onclick) #Signal

        #Button to start the alignement
        self.btnSACalculate = QPushButton("Align")
        tabSliceAlignment_l.addWidget(self.btnSACalculate)
        self.btnSACalculate.clicked.connect(self.btnSACalculate_onclick) #Signal


    def ui_do_tab_ResMeas(self):
        # /tabQuoll \
        if quoll.is_quoll_available:
        #if True: #Debug UI in windows
            tabQuoll = QWidget()
            self.tabwidget.addTab(tabQuoll, "Res. measurement")

            vBox4 = QVBoxLayout() #Items organised vertically
            tabQuoll.setLayout(vBox4)

            label0 = QLabel("Resolution measurement\nFourier Ring Correlation (FRC) with Quoll")
            label0.setWordWrap(True)
            vBox4.addWidget(label0)

            qgrid3_layout = QGridLayout()
            vBox4.addLayout(qgrid3_layout)
            qgrid3_layout.addWidget(QLabel("tile size (px):"),0,0)
            self.spbQuolTileSize_px = QSpinBox()
            self.spbQuolTileSize_px.setMinimum(128) #miplib minimum is 128 pixels
            self.spbQuolTileSize_px.setMaximum(10000)
            self.spbQuolTileSize_px.setValue(256)
            qgrid3_layout.addWidget(self.spbQuolTileSize_px, 0,1)
            
            #add pixel size box
            qgrid3_layout.addWidget(QLabel("pixel size (nm):"),1,0)
            self.spbQuolPixelSize_nm = QDoubleSpinBox()
            self.spbQuolPixelSize_nm.setMinimum(1e-6)
            self.spbQuolPixelSize_nm.setMaximum(10000000)
            self.spbQuolPixelSize_nm.setValue(1.0)
            qgrid3_layout.addWidget(self.spbQuolPixelSize_nm, 1,1)           

            btnQuollCalcFRCTiled = QPushButton("Calculate FRC tiled")
            vBox4.addWidget(btnQuollCalcFRCTiled)
            btnQuollCalcFRCTiled.clicked.connect(self.btnQuollCalcFRCTiled_onclick) #Signal

            # add option to choose calibration function
            qgrid3_layout.addWidget(QLabel("calibration function:"), 2, 0)
            from inspect import getmembers, isfunction
            self.spbQuollCF = QComboBox()
            self.spbQuollCF.addItem("None") # do not apply calibration
            from quoll.frc import frc_calibration_functions
            cfs = getmembers(frc_calibration_functions, isfunction)
            self.spbQuollCF.addItems([cfs[i][0] for i in range(len(cfs))])
            qgrid3_layout.addWidget(self.spbQuollCF, 2, 1)

            #table with infomation
            quollTableContainerWidget=QWidget()
            vBox4.addWidget(quollTableContainerWidget)
            quollTableContainerWidgetLayout=QVBoxLayout()
            quollTableContainerWidget.setLayout(quollTableContainerWidgetLayout)
            self.quollTable=widget_OkapiemTable(self.viewer)
            quollTableContainerWidgetLayout.addWidget(self.quollTable)

    def ui_do_tab_Filters(self):
        tabFilters = QWidget()
        self.tabwidget.addTab(tabFilters, "Filters")

        vboxlay = QVBoxLayout() #Items organised vertically
        tabFilters.setLayout(vboxlay)

        label0 = QLabel("Filters")
        vboxlay.addWidget(label0)

        #Optical flow denoise start
        self.rbFlowDenoise= QRadioButton("Optical Flow Denoising")
        self.rbFlowDenoise.setChecked(True)
        vboxlay.addWidget(self.rbFlowDenoise)

        grpBox1 = QGroupBox("") #Box, no title, radiobutton does it
        vboxlay.addWidget(grpBox1)
        vbox2 = QVBoxLayout()
        grpBox1.setLayout(vbox2)
        label1 = QLabel("Parameters")
        vbox2.addWidget(label1)

        #Grid for sigma z,y,x
        qgrid2_layout = QGridLayout()
        #qgrid2_layout.setAlignment(Qt.AlignmentFlag.AlignRight) # does not work
        vbox2.addLayout(qgrid2_layout)
        qgrid2_layout.addWidget(QLabel(""),0,0)
        qgrid2_layout.addWidget(QLabel("z"),0,1, Qt.AlignmentFlag.AlignHCenter)
        qgrid2_layout.addWidget(QLabel("x"),0,3, Qt.AlignmentFlag.AlignHCenter)
        qgrid2_layout.addWidget(QLabel("y"),0,2, Qt.AlignmentFlag.AlignHCenter)
        qgrid2_layout.addWidget(QLabel("sigma"),1,0, Qt.AlignmentFlag.AlignRight)
        self.dspb_OF_sigma_z=QDoubleSpinBox()
        self.dspb_OF_sigma_z.setMinimum(0.1)
        self.dspb_OF_sigma_z.setMaximum(1000.0)
        self.dspb_OF_sigma_z.setDecimals(2)
        self.dspb_OF_sigma_z.setSingleStep(0.1)
        self.dspb_OF_sigma_z.setValue(2.0)
        qgrid2_layout.addWidget(self.dspb_OF_sigma_z,1,1)

        self.dspb_OF_sigma_y=QDoubleSpinBox()
        self.dspb_OF_sigma_y.setMinimum(0.1)
        self.dspb_OF_sigma_y.setMaximum(1000.0)
        self.dspb_OF_sigma_y.setDecimals(2)
        self.dspb_OF_sigma_y.setSingleStep(0.1)
        self.dspb_OF_sigma_y.setValue(2.0)
        qgrid2_layout.addWidget(self.dspb_OF_sigma_y,1,2)

        self.dspb_OF_sigma_x=QDoubleSpinBox()
        self.dspb_OF_sigma_x.setMinimum(0.1)
        self.dspb_OF_sigma_x.setMaximum(1000.0)
        self.dspb_OF_sigma_x.setDecimals(2)
        self.dspb_OF_sigma_x.setSingleStep(0.1)
        self.dspb_OF_sigma_x.setValue(2.0)
        qgrid2_layout.addWidget(self.dspb_OF_sigma_x,1,3)

        qgrid2_layout.addWidget(QLabel("levels"),2,0, Qt.AlignmentFlag.AlignRight)
        self.spb_OF_levels=QSpinBox()
        self.spb_OF_levels.setMinimum(0)
        self.spb_OF_levels.setMaximum(100)
        self.spb_OF_levels.setValue(3)
        self.spb_OF_levels.setSingleStep(1)
        qgrid2_layout.addWidget(self.spb_OF_levels,2,1)

        qgrid2_layout.addWidget(QLabel("winsize"),2,2, Qt.AlignmentFlag.AlignRight)
        self.spb_OF_winsize=QSpinBox()
        self.spb_OF_winsize.setMinimum(0)
        self.spb_OF_winsize.setMaximum(100)
        self.spb_OF_winsize.setValue(5)
        self.spb_OF_winsize.setSingleStep(1)
        qgrid2_layout.addWidget(self.spb_OF_winsize,2,3)

        qgrid2_layout.addWidget(QLabel("p-processes"),3,0, Qt.AlignmentFlag.AlignRight)
        self.spb_OF_maxproc=QSpinBox()
        self.spb_OF_maxproc.setMinimum(0)
        self.spb_OF_maxproc.setMaximum(100)
        self.spb_OF_maxproc.setValue(4)
        self.spb_OF_maxproc.setSingleStep(1)
        qgrid2_layout.addWidget(self.spb_OF_maxproc,3,1)

        qgrid2_layout.addWidget(QLabel("mode"),3,2, Qt.AlignmentFlag.AlignRight) 
        self.cb_OF_procmode = QComboBox()
        self.cb_OF_procmode.addItem("threaded")
        self.cb_OF_procmode.addItem("multiproc")
        self.cb_OF_procmode.addItem("sequential")
        qgrid2_layout.addWidget(self.cb_OF_procmode,3,3)        

        #Button to start the alignement
        self.btn_Filters_Filter = QPushButton("Filter")
        vboxlay.addWidget(self.btn_Filters_Filter)
        self.btn_Filters_Filter.clicked.connect(self.btn_Filters_Filter_onclick) #Signal
        #Optical flow denoise end

    def freqselector_high_valueChanged(self):
        #check value is not lower than lowvalue
        v0 = self.freqselector_high.value()
        if v0< self.freqselector_low.value():
            self.freqselector_high.setValue( 0)
        #self.doFFTDirFilterSlice()
    
    def freqselector_low_valueChanged(self):
        v0 = self.freqselector_low.value()
        if v0> self.freqselector_high.value():
            self.freqselector_low.setValue(v0-1)
        #self.doFFTDirFilterSlice()

    def doFFTDirFilterSlice(self):
        
        print("doFFTDirFilterSlice()")

        curDataSlice = np.asarray(self.w_setselect.get_active_image_selected_data_slice())
        if not curDataSlice is None:

            datares=None

            #Check which filter to use
            if self.rbDirFFT.isChecked():
                freq_low = self.freqselector_low.value()
                freq_high = self.freqselector_high.value()
                
                datares = filters.fft_bandpass_filter_dirx_2D(curDataSlice,freq_low, freq_high)

            elif self.rbChafer.isChecked():
                curLabelsSlice= np.asarray(self.w_setselectlabel.get_active_image_selected_data_slice())
                c_filter = chafer.cls_charge_artifact_suppression_filter(
                    self.chafer_widg_nlines.value(),
                    self.chafer_widg_maxlength.value(),
                    self.chafer_widg_minlength.value()
                    )
                
                datares, opts = c_filter.charge_artifact_FD_filter_downup_av_prevlines3_2d(curDataSlice, curLabelsSlice)

            if not datares is None:
                #Not sure if using a single preview data layer is ok
                if self.myFFTFilterDirSliceRes is None or not self.myFFTFilterDirSliceRes in self.viewer.layers:
                    self.myFFTFilterDirSliceRes = self.viewer.add_image(datares, name="c-filtered", contrast_limits= (curDataSlice.min(), curDataSlice.max()))
                else:
                    #Dont create a new layer, just modify it
                    self.myFFTFilterDirSliceRes.data = datares
                    #self.myImageLayer.refresh()

    def doFFTDirFilterWhole(self):
        #From current layer/data selected, apply the filter to the whole data
        
        print("doFFTDirFilterWhole()")
        curData = np.asarray(self.w_setselect.get_active_image_selected_data())

        if not curData is None:
            
            data_filt_whole=None
            
            cCalculate=None

            #Check which filter to use
            if self.rbDirFFT.isChecked():
                freq_low = self.freqselector_low.value()
                freq_high = self.freqselector_high.value()
            
                def func(data2d):
                    return filters.fft_bandpass_filter_dirx_2D(data2d,freq_low, freq_high)

                cCalculate = cSliceBySliceFunctHelper(self.viewer,func, curData)

            
            elif self.rbChafer.isChecked():
                #Use chafer to remove charging artifacts
                #Use charge artifact label

                #Dont use the Chafer's 3d method
                #Use this engine to do slice by slice so that it can be monitored

                curLabels= np.asarray(self.w_setselectlabel.get_active_image_selected_data())

                c_filter = chafer.cls_charge_artifact_suppression_filter(
                    self.chafer_widg_nlines.value(),
                    self.chafer_widg_maxlength.value(),
                    self.chafer_widg_minlength.value()
                    )

                def func(data2d, labels2d):
                    datares, opts = c_filter.charge_artifact_FD_filter_downup_av_prevlines3_2d(data2d, labels2d)
                    return datares

                cCalculate = cSliceBySliceFunctHelper(self.viewer,func, curData, curLabels)
            
            if not cCalculate is None:

                estimate_iters = cCalculate.get_len()
                print(f"estimate_iters:{estimate_iters}")
                self.iter=0

                pbr=progress(total=estimate_iters,desc="Chafer charging artifact removal in progress")

                def calbackfn():
                    print(f"slice iter:{self.iter} / {estimate_iters}")
                    self.iter+=1
                    pbr.update(1)
                    pbr.refresh()

                data_filt_whole = cCalculate.run(calbackfn)

                print("Calculation complete")

                pbr.close()

                if not data_filt_whole is None:
                    self.viewer.add_image(data_filt_whole, name="chafer filter")
    
    def btnMCCalculate_onclick(self):
        #Run the measure charging

        print("btnMCCalculate_onclick()")

        #Read tile size parameter and gets current image
        tilesize = self.mc_tilesize_selector.value()

        data0 = np.asarray(self.w_setselect.get_active_image_selected_data_slice())

        curDataSlice=data0
        if not curDataSlice is None:
            #Calculate
            overlay = measure_charging.generate_heatmap(curDataSlice, tilesize)
            self.viewer.add_image(overlay, name="Charging artefacts", opacity=0.3, colormap="viridis")


    def btnSACalculate_onclick(self):
        print("btnSACalculate_onclick()")

        res=None

        data3d = self.w_setselect.get_active_image_selected_data()
        if  data3d is None:
            raise ValueError("Could not get active image data.")

        if data3d.ndim!=3:
            raise ValueError("Data selected is not 3D")
        

        #Ensure data3d is napari and not dask
        data3d = np.asarray(data3d)

        #Estimate number of iterations for progress bar
        nslices = data3d.shape[0]
        niterations = 2*(nslices+1) #TODO: Proably should be 2nslices

        pbr=progress(total=niterations,desc="Alignment progress")

        def callbkfn():
            pbr.update(1)
            pbr.refresh()

        #Show or not show Activity dialog?
        #activ_dialog=self.viewer.window._qt_viewer.window()._activity_dialog #Warning this will be unavailable in the future
        #activ_dialog.show()

        #res= slice_alignment.align_stack(data3d, slice_alignment.ALIGNMENT_METHOD_DEFAULT,callbkfn)
        # sa_method = slice_alignment.ALIGNMENT_METHOD_DEFAULT

        # sa_method['translation'] = self.chkbxTranslate.isChecked()
        # #sa_method['affine'] = self.chkbxAffine.isChecked()
        # sa_method['linear'] = self.chkbxLinear.isChecked()
        # sa_method['rotation'] = self.chkbxRotation.isChecked()
        # sa_method['shearing_x'] = self.chkbxShearX.isChecked()
        # sa_method['shearing_y'] = self.chkbxShearY.isChecked()
        # sa_method['scaling'] = self.chkbxScaling.isChecked()
        # sa_method['stretching_x'] = self.chkbxStretchX.isChecked()
        # sa_method['stretching_y'] = self.chkbxStretchY.isChecked()

        #res= slice_alignment.align_stack(data3d, sa_method,callbkfn)

        p0= self.get_UI_SA_parameters()

        res, _= slice_alignment.align_stack1(data3d,
                                        p0['method'],
                                        d2dmedian_ratio=p0['d2dmedian'],
                                        ratio_test=p0['reldist2match'],
                                        callbk_tick_fn=callbkfn
                                    )

        pbr.close()
        #activ_dialog.hide()

        if not res is None:
            self.viewer.add_image(res, name="stack aligned")
        

    def btnQuollCalcFRC_onclick(self):
        #NOT USED

        print("btnQuollCalcFRC_onclick()")

        with progress(total=2,desc="FRC calculation in progress") as pbr:

            data2d = np.asarray(self.w_setselect.get_active_image_selected_data_slice())

            pbr.update(1)

            res = quoll.getFRC(data2d)

            pbr.update(1)

            print (f"FRC resolution : {res}")

            #print result on a widget
            self.lblFRCWholeImg.setText( f"{res:.4e}")
        return


    def btnQuollCalcFRCTiled_onclick(self):
        print("btnQuollCalcFRCTiled_onclick()")

        #Gathers tile size parameter
        tile_size = self.spbQuolTileSize_px.value()
        pixel_size_nm = self.spbQuolPixelSize_nm.value()
        cf_in = str(self.spbQuollCF.currentText())
        
        res = None

        #Runs the FRC calculation for tiles
        with progress(total=2,desc="FRC tiled calculation in progress") as pbr:

            data2d = np.asarray(self.w_setselect.get_active_image_selected_data_slice())

            pbr.update(1)

            res = quoll.getTiledFRC(
                data2d, 
                cf_in=cf_in,
                tilesize=tile_size, 
                pixel_size_nm=pixel_size_nm
            )

            pbr.update(1)

            print (f"FRC tiled result : {res}")

            if not res is None:
                dfres, heatmap = res
                print(f"FRC result summary:{dfres.describe().to_dict()}")
                self.quollTable.setFromDict(dfres.describe().to_dict()['Resolution'])

                # display a heatmap
                self.viewer.add_image(heatmap,name="FRC heatmap",
                    colormap="viridis",
                    opacity=0.3
                )
        return

    def btn_testSA_onclick(self):
        #Preview of slice alignment with current slice and next

        #Code similar to get_active_image_selected_data_slice()
        curData = self.w_setselect.get_active_image_selected_data()

        if curData is None:
            raise ValueError("No data")
        
        if curData.ndim!=3:
            raise ValueError("Data is not 3-dimensional. This plugin only works with 3D data")
        
        dims0 = self.viewer.dims
        nslice = dims0.current_step[0] #Gets the slice number that is currently being viewed, 1st index is the slice number
        if nslice==curData.shape[0]:
            raise ValueError("Cannot use this button with the last slice.")
        
        order0 = dims0.order[0] #Gets the axis (or plane), that the Viewer is currently viewing
        if order0 != 0: #z
            ValueError("Stack alignment only works along z direction. Please change visualisation to XY plane and ndisplay=2 for this button to work.")

        data_2slice = np.asarray(curData[nslice:nslice+2, :,:])

        p0= self.get_UI_SA_parameters()

        res, _= slice_alignment.align_stack1(data_2slice,
                                          p0['method'],
                                          d2dmedian_ratio=p0['d2dmedian'],
                                          ratio_test=p0['reldist2match']
                                          )

        if not res is None:
            self.viewer.add_image(res, name="stack-align preview")


    def btn_Filters_Filter_onclick(self):
        #Read parameters and run OF
        
        if self.rbFlowDenoise.isChecked():
            #Flow denoise filtering

            with progress(range(3)) as pbr:
                try:
                    import flowdenoising.flowdenoising_mod as fdn
                except ModuleNotFoundError:
                    
                    raise ModuleNotFoundError(
                            """flowdenoising (rfi fork) package not installed.
                            Please install it by running the following command in the napari console:
                            > pip install git+https://github.com/rosalindfranklininstitute/FlowDenoising.git"""
                        )

                sigma_z = self.dspb_OF_sigma_z.value()
                sigma_y = self.dspb_OF_sigma_y.value()
                sigma_x = self.dspb_OF_sigma_x.value()
                sigma0 = (sigma_z, sigma_y, sigma_x)

                levels0=self.spb_OF_levels.value()
                winsize0=self.spb_OF_winsize.value()
                nprocesses0=self.spb_OF_maxproc.value()

                process_mode0 = self.cb_OF_procmode.currentText()
                
                data0=np.asarray(self.w_setselect.get_active_image_selected_data()) #get active data


                #Run the OF flow here
                of_filter = fdn.cFlowDenoiser(
                    sigma=sigma0,
                    levels=levels0,
                    winsize=winsize0,
                    max_number_of_processes=nprocesses0,
                    process_mode=process_mode0)
                
                
                pbr.update(1)
                
                res0 = of_filter.runOpticalFlow(data0)

                pbr.update(1)

                #Add filtered image
                if not res0 is None:
                    self.viewer.add_image(res0,
                                        name="flow denoised")
                

    def get_UI_SA_parameters(self):
        sa_method = slice_alignment.ALIGNMENT_METHOD_DEFAULT
        sa_method['translation'] = self.chkbxTranslate.isChecked()
        #sa_method['affine'] = self.chkbxAffine.isChecked()
        sa_method['linear'] = self.chkbxLinear.isChecked()
        sa_method['rotation'] = self.chkbxRotation.isChecked()
        sa_method['shearing_x'] = self.chkbxShearX.isChecked()
        sa_method['shearing_y'] = self.chkbxShearY.isChecked()
        sa_method['scaling'] = self.chkbxScaling.isChecked()
        sa_method['stretching_x'] = self.chkbxStretchX.isChecked()
        sa_method['stretching_y'] = self.chkbxStretchY.isChecked()

        d2dmedian = self.spb_SliceAlign_D2Dmedian_Ratio.value()
        reldist2match = self.spb_SliceAlign_RatioTest.value()

        return {'method':sa_method,'d2dmedian':d2dmedian, 'reldist2match':reldist2match}
                

    def hello_World(self):
        #prints on the console something enlighting. Only used in CI automated tests
        print ("Hello World, from Okapi-EM!")


class widget_ImageSetCurrentSelect(QWidget):
    '''
    A QtWidget that helps select Image data for further processing

    Needs to be initialised with a napari Viewer object

    '''

    def __init__(self, napari_viewer: 'napari.viewer.Viewer', header=""):
        super().__init__()

        self.viewer = napari_viewer #Reference to the viewer that will be needed later
        #self.viewer = napari.current_viewer() #Doesnt work, returns None

        layout0 = QVBoxLayout()
        self.setLayout(layout0)

        group_box0 = QGroupBox(header)
        #self.t0_vbox0.addWidget(self.group_box0)
        layout0.addWidget(group_box0)

        gridlayout = QGridLayout()
        group_box0.setLayout(gridlayout)

        self.b_set = QPushButton("Set currently selected")
        gridlayout.addWidget(self.b_set,0,0)
        self.b_set.clicked.connect(self.b_set_on_click) #Signal
  
        gridlayout.addWidget(QLabel("name"),1,0)
        self.l_name_v = QLabel("")
        gridlayout.addWidget(self.l_name_v,1,1)

        gridlayout.addWidget(QLabel("shape"),2,0)
        self.l_shape_v = QLabel("")
        gridlayout.addWidget(self.l_shape_v,2,1)

        #layout.addWidget(QLabel("slicing axis"),3,0)
        self.l_slicingaxis_v = QLabel("")
        # layout.addWidget(self.l_slicingaxis_v,3,1)
        # #self.l_slicingaxis_v.hide() #Hide
        # self.l_slicingaxis_v.setHidden(True) #Hide

        #layout.addWidget(QLabel("slice shape"),4,0)
        self.l_sliceshape_v = QLabel("")
        # layout.addWidget(self.l_sliceshape_v,4,1)
        # self.l_sliceshape_v.hide() #Hide

        #layout.addWidget(QLabel("nslice"),5,0)
        self.l_nslice_v = QLabel("")
        # layout.addWidget(self.l_nslice_v,5,1)
        # self.l_nslice_v.hide() #Hide

        self.curImage=None

    def b_set_on_click(self):
        print("b_set_on_click()")
        self.set_active_image_from_nap_layers()
        
        #active0=self.get_active_image() #This also sets default image   
        return

    def set_default_napari_image(self,nap_image):
        print("set_default_napari_image()")
        self.curImage=None
        print(f"type(nap_image) : {type(nap_image) }")
        if isinstance(nap_image, napari.layers.Image) or isinstance(nap_image, napari.layers.Labels):
            self.l_name_v.setText(nap_image.name)
            print("Current active layer name is: ", nap_image.name)
            self.curImage = nap_image
            self.l_shape_v.setText(str(self.curImage.data.shape))
        else:
            print("Not a valid layer for Okapi-EM.")
            self.curImage=None
            self.l_name_v.setText("")
            self.l_shape_v.setText("")
        
        return

    def get_active_image(self):
        print("get_active_image()`")
        
        if not self.curImage is None:
            if self.curImage in self.viewer.layers: #Ensure that image was not deleted
                return self.curImage
                
        #curImage not available anymore

        # try to grab an image from selected
        self.set_active_image_from_nap_layers()
        if not self.curImage is None:
            return self.curImage

        return None
    
    def set_active_image_from_nap_layers(self):
        #Set active image from napari layers
        print("set_active_image_from_nap_layers()")
        
        active0=None
        if not self.viewer is None:
            active0 = self.viewer.layers.selection.active #check if any layer is active
            self.set_default_napari_image(active0)

        else:
            raise ValueError("No Viewer")

        return

    def get_active_image_selected_data(self) ->ImageData:
        print("get_active_image_selected_data()")

        im = self.get_active_image()
        if not im is None:
            return self.curImage.data
        
        return None
    
        # if self.curImage in self.viewer.layers: #Ensure that image was not deleted
        #     return self.curImage.data
        # else:
            
        
            # if self.curImage is None:
            #     #Try to grab an image
            #     self.b_set_on_click() #Does the same thing as clicking the button to set default image data
                
            #     if self.curImage is None: #Check if any layer is selected
            #         #Clear selection

            #         raise ValueError("Could not get active image data.")
            #     else:
                    


    def get_active_image_selected_data_slice(self):
        print("get_active_image_selected_data_slice()")

        curDataSlice=None

        curData = self.get_active_image_selected_data()
        #print(f"curData:{curData}")

        if not curData is None:
            if curData.ndim==3:
                if self.viewer.dims.ndisplay == 2:
                    #Get the currently viewing slice
                    dims0 = self.viewer.dims

                    nslice = dims0.current_step[0] #Gets the slice number that is currently being viewed, 1st index is the slice number, but which axis is not given here
                    print("nslice: ", str(nslice))
                    self.l_nslice_v.setText(str(nslice))

                    order0 = dims0.order[0] #Gets the axis (or plane), that the Viewer is currently viewing
                    if order0 == 0: #z
                        curDataSlice = self.curImage.data[nslice,...]
                        self.l_slicingaxis_v.setText('z')
                    elif order0 == 1: #y
                        curDataSlice = self.curImage.data[:,nslice,:]
                        self.l_slicingaxis_v.setText("y")
                    else: #x
                        curDataSlice = self.curImage.data[...,nslice]
                        self.l_slicingaxis_v.setText("x")

                    #Slice shape
                    s_shape = curDataSlice.shape
                    self.l_sliceshape_v.setText(str(s_shape))

                else:
                    print("Viewer is in 3D mode, cannot get slice")
                    self.l_slicingaxis_v.setText("Viewer in 3D mode, cannot get slice")
            else:
                print("ndim=2. The image layer itself is the slice.")
                curDataSlice= curData
            
        return curDataSlice
        

class cSliceBySliceFunctHelper():

    def __init__(self,
            napari_viewer: 'napari.viewer.Viewer',
            func,
            data, 
            data_opt=None 
        ):
        '''
            data
            func to apply to slices. Must be like func(data2D) -> datares2D

            Some functions require two datasets, hence the optional data3d_opt
            In this case func(data2D, data3d_opt) -> datares2D
        '''

        #viewer = napari.current_viewer()
        viewer = napari_viewer

        self.datadim = data.ndim
        self.data_oriented = None
        self.data_opt_oriented = None

        if data.ndim==3:
            #gets orientation of viewing. If viewer is in 3D it sets to z-axis order=0
            self.order=0
            if viewer.dims.ndisplay == 2:
                dims0 = viewer.dims
                self.order = dims0.order[0] #Gets the axis (or plane), that the Viewer is currently viewing

            self.has_opt = False
            if not data_opt is None:
                self.has_opt = True

            if self.order == 0: #z
                self.data_oriented=data
                #self.l_slicingaxis_v.setText('z')
                if self.has_opt:
                    self.data_opt_oriented=data_opt

            elif self.order == 1: #y
                self.data_oriented = np.transpose(data,(1,2,0))
                #self.l_slicingaxis_v.setText("y")
                if self.has_opt:
                    self.data_opt_oriented= np.transpose(data_opt,(1,2,0))
            else: #x
                self.data_oriented = np.transpose(data,(2,0,1))
                #self.l_slicingaxis_v.setText("x")
                if self.has_opt:
                    self.data_opt_oriented= np.transpose(data_opt,(2,0,1))

            self.n= self.data_oriented.shape[0]

            self.func = func

        else:
            print("data is not 3D. Function will be applied to the single slice")
            self.data_oriented = data
            self.data_opt_oriented = data_opt
            self.n=1
    
    def run(self, callbackfn=None):
        if self.datadim==3:
            datares_shaperestored=None
            
            if not self.data_oriented is None:
                
                #Check result type with first slice
                if not self.has_opt:
                    slice0res = self.func(self.data_oriented[0,:,:])
                else:
                    slice0res = self.func(self.data_oriented[0,:,:], self.data_opt_oriented[0,:,:])

                #Set result to have the same type as first layer result
                datares = np.zeros_like(self.data_oriented, dtype=slice0res.dtype)

                datares[0,:,:]= slice0res

                #callback
                if not callbackfn is None:
                    callbackfn()

                #Do for the rest of slices
                for i in range(1, self.n):
                    #dataslice = self.data_oriented[i,:,:]

                    #Apply function
                    if not self.has_opt:
                        datasliceres = self.func(self.data_oriented[i,:,:])
                    else:
                        datasliceres = self.func(self.data_oriented[i,:,:],self.data_opt_oriented[i,:,:])
                    
                    datares[i,:,:] = datasliceres

                    #callback
                    if not callbackfn is None:
                        callbackfn()
                    
                #Upon completion restore the shape of the data
                if self.order == 0: #z
                    datares_shaperestored = datares
                elif self.order == 1: #y
                    datares_shaperestored = np.transpose(datares,(2,0,1))
                else: #x
                    datares_shaperestored = np.transpose(datares,(1,2,0))

            return datares_shaperestored
        
        elif self.datadim==2:
            data2dres = None
             #Apply function
            if not self.has_opt:
                data2dres = self.func(self.data_oriented)
            else:
                data2dres = self.func(self.data_oriented,self.data_opt_oriented)
            return data2dres
        
        return None
    
    def get_len(self):
        #Returns number of iterations expected for the calculation

        niter=0
        if not self.data_oriented is None:
            niter= self.n

        return niter


class widget_OkapiemTable(QTableWidget):
    '''
    A simple table that displays a pandas dataframe without interaction

    '''
    def __init__(self, napari_viewer: 'napari.viewer.Viewer'):
        super().__init__()

    def setDataframe(self,dataframe):
        self.clearContents()

        self.dataframe = dataframe

        #Fill table with content from dataframe
        rowcount = dataframe.shape[0]+1
        columncount = dataframe.shape[1]

        self.setRowCount(rowcount)
        self.setColumnCount(columncount)

        #Header
        for i,colnames in enumerate(dataframe.keys()):
            self.setItem(0,i, QTableWidgetItem(str(colnames)))
        
        #data
        for i, colnames in enumerate(dataframe.keys()):
            for j, value in enumerate(dataframe.get(colnames)):
                self.setItem(j, i, QTableWidgetItem(str(value)))

    def setFromDict(self, dict0,numberformat='{:.3e}'):
        self.clearContents()

        rowcount = len(dict0)
        columncount= 2
        self.setRowCount(rowcount)
        self.setColumnCount(columncount)
        
        self.horizontalHeader().setVisible(False)
        self.verticalHeader().setVisible(False)

        for i, key0 in enumerate(dict0.keys()):
            v= dict0[key0]
            self.setItem(i, 0, QTableWidgetItem(str(key0)))
            #self.setItem(i, 1, QTableWidgetItem(str(v)))
            self.setItem(i, 1, QTableWidgetItem(numberformat.format(v)))