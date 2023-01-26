#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pickle
import sys
import os
import os.path

#scale QT applications for high resolution screens 
qtscale = float(os.environ.get('QT_SCALE_FACTOR',1))


#!/usr/bin/env python
# -*- coding: utf-8 -*-
try:
    from PyQt5.QtCore import *
    from PyQt5.QtGui import *
    from PyQt5.QtWidgets import *
except:
    from PyQt4.QtCore import *
    from PyQt4.QtGui import *


sys.path.append('../')
from make_graphs import *
from sawtooths import sawtooths_detection
from imp_analysis import imp_analysis
import traceback
from prepare_data import *
from main import tomography
import pytomo
import time
from shared_modules import debug
import gc
import config
from matplotlib import rcParams

"""
Graphical user interface for Tomography
subprograms: make_graphs - create output graphs and graphs for GUI
             prepare_data - loads data of from file/internet, prepare them to correct format
             start.py    -  main core program, calls all other
inputs:             tomography.npy - file with saved settings, will be automaticly removed when coputer is changed
Tomáš Odstrčil,(Michal Odstrčil), IPP, 2015
"""



class MainWindow(QMainWindow):
    def __init__(self, setting, tokamak):
        QMainWindow.__init__(self, None)
        


        self.setting = setting
        self.tokamak = tokamak

        self.centralwidget = QWidget(self)
        self.gridLayout = QGridLayout(self.centralwidget)
        self.setWindowTitle('Tomography')

        self.Expand = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.Fixed  = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
 

        self.main_tab = QTabWidget(self)
        self.main_tab.setTabPosition(QTabWidget.West)


        self.gridLayout.addWidget(self.main_tab, 0, 0, 1,2)

        self.tables_dict = {}
        
        # IMAGES IN TABLE !!!
        i_tab = 0
        self.tab_widget_reconstruction = QWidget(self)
        self.main_tab.insertTab(i_tab,self.tab_widget_reconstruction, 'Main')
        self.tables_dict['Main'] = i_tab
        i_tab+=1
        self.tab_widget_reconstruction.setSizePolicy(self.Expand)
        self.verticalLayout_rec = QVBoxLayout(self.tab_widget_reconstruction)
        self.horizontalLayout_rec = QHBoxLayout()
        #preview right
        self.Picture_recons_1 = QLabel()
        #self.Picture_recons_1.setText('Preview profile X')
        self.Picture_recons_1.setAlignment(Qt.AlignCenter)
        self.Picture_recons_1.setScaledContents(True)
        self.Picture_recons_1.setSizePolicy(self.Fixed)
        self.horizontalLayout_rec.addWidget(self.Picture_recons_1)


        #preview left
        self.Picture_recons_2 = QLabel()
        self.Picture_recons_2.setAlignment(Qt.AlignCenter)
        self.Picture_recons_2.setScaledContents(True)
        self.Picture_recons_2.setSizePolicy(self.Fixed)
        self.horizontalLayout_rec.addWidget(self.Picture_recons_2)
        self.verticalLayout_rec.addLayout(self.horizontalLayout_rec)
        # slider
        self.slider_rec = QSlider(Qt.Horizontal)
        self.slider_rec_spin = QSpinBox()
        self.slider_rec.setVisible(False)
        self.slider_rec_spin.setVisible(False)

        self.horizontalLayout_9 = QHBoxLayout()
        self.horizontalLayout_9.addWidget(self.slider_rec)
        self.horizontalLayout_9.addWidget(self.slider_rec_spin)
        self.verticalLayout_rec.addLayout(self.horizontalLayout_9)



        # overview plot
        #TODO only for AUG now!!!
        if os.path.exists('geometry/'+tokamak.name+'/overview_plots/'):
            self.tab_widget_overview = QWidget(self)
            self.main_tab.insertTab(i_tab,self.tab_widget_overview, 'Overview')
            self.main_tab.setTabEnabled(i_tab,True)
            self.tables_dict['Overview'] = i_tab
            i_tab+=1
            self.OverviewPlot = QLabel()
            self.OverviewPlot.setAlignment(Qt.AlignCenter)
            self.change_overview_plots()
            self.horizontalLayout_over = QHBoxLayout()
            self.verticalLayout_over = QVBoxLayout(self.tab_widget_overview)
            self.horizontalLayout_over.addWidget(self.OverviewPlot)
            self.verticalLayout_over.addLayout(self.horizontalLayout_over)


        # SVD
        self.tab_widget_SVD = QWidget(self)
        self.verticalLayout_SVD = QVBoxLayout(self.tab_widget_SVD)
        self.horizontalLayout_SVD = QHBoxLayout()
        self.main_tab.insertTab(i_tab,self.tab_widget_SVD, 'SVD analysis')
        self.main_tab.setTabEnabled(i_tab,False)
        self.tables_dict['SVD'] = i_tab
        i_tab+=1
        
        #SVD right
        self.Picture_SVD_1 = QLabel()
        self.Picture_SVD_1.setAlignment(Qt.AlignCenter)
        self.Picture_SVD_1.setScaledContents(True)
        self.Picture_SVD_1.setSizePolicy(self.Fixed)
        self.horizontalLayout_SVD.addWidget(self.Picture_SVD_1)
        #SVD left

        self.Picture_SVD_2 = QLabel()
        self.Picture_SVD_2.setAlignment(Qt.AlignCenter)
        self.Picture_SVD_2.setScaledContents(True)
        self.Picture_SVD_2.setSizePolicy(self.Fixed)
        self.horizontalLayout_SVD.addWidget(self.Picture_SVD_2)
        self.verticalLayout_SVD.addLayout(self.horizontalLayout_SVD)
        

        # slider
        self.slider_SVD = QSlider(Qt.Horizontal)
        self.slider_SVD_spin = QSpinBox()
        self.slider_SVD.setVisible(False)
        self.slider_SVD_spin.setVisible(False)
        self.horizontalLayout_SVD = QHBoxLayout()
        self.horizontalLayout_SVD.addWidget(self.slider_SVD)
        self.horizontalLayout_SVD.addWidget(self.slider_SVD_spin)
        self.verticalLayout_SVD.addLayout(self.horizontalLayout_SVD)



        # position
        self.tab_widget_position= QWidget(self)
        self.gridLayout_position= QGridLayout(self.tab_widget_position)
        self.main_tab.insertTab(i_tab,self.tab_widget_position, 'Position')
        self.main_tab.setTabEnabled(i_tab,False)
        self.tables_dict['Position'] = i_tab
        i_tab+=1

        #postprocess right
        self.Picture_position_1 = QLabel()
        self.Picture_position_1.setAlignment(Qt.AlignCenter)
        self.Picture_position_1.setScaledContents(True)
        self.Picture_position_1.setSizePolicy(self.Fixed)
        self.gridLayout_position.addWidget(self.Picture_position_1, 0, 0)

        #postprocess left
        self.Picture_position_2 = QLabel()
        self.Picture_position_2.setAlignment(Qt.AlignCenter)
        self.Picture_position_2.setScaledContents(True)
        self.Picture_position_2.setSizePolicy(self.Fixed)
        self.gridLayout_position.addWidget(self.Picture_position_2, 0, 1)




        # emissivity
        self.tab_widget_emis= QWidget(self)
        self.gridLayout_emis= QGridLayout(self.tab_widget_emis)
        self.main_tab.insertTab(i_tab,self.tab_widget_emis, 'Emissivity')
        self.main_tab.setTabEnabled(i_tab,False)
        self.tables_dict['Emissivity'] = i_tab
        i_tab+=1


        # left
        self.Picture_emis_1 = QLabel()
        self.Picture_emis_1.setAlignment(Qt.AlignCenter)
        self.Picture_emis_1.setScaledContents(True)
        sizePolicy = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        self.Picture_emis_1.setSizePolicy(sizePolicy)
        self.gridLayout_emis.addWidget(self.Picture_emis_1, 0, 0)

        # right
        self.Picture_emis_2 = QLabel()
        self.Picture_emis_2.setAlignment(Qt.AlignCenter)
        self.Picture_emis_2.setScaledContents(True)
        sizePolicy = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        self.Picture_emis_2.setSizePolicy(sizePolicy)
        self.gridLayout_emis.addWidget(self.Picture_emis_2, 0, 1)




        # advanced
        if setting['post_proces_equi']:
            self.tab_widget_advanced= QWidget(self)
            self.gridLayout_advanced= QGridLayout(self.tab_widget_advanced)
            self.main_tab.insertTab(i_tab,self.tab_widget_advanced, 'Equilibrium')
            self.main_tab.setTabEnabled(i_tab,False)
            self.tables_dict['Equilibrium'] = i_tab
            i_tab+=1


            #advanced left
            self.Picture_advanced_1 = QLabel()
            self.Picture_advanced_1.setAlignment(Qt.AlignCenter)
            self.Picture_advanced_1.setScaledContents(True)
            self.Picture_advanced_1.setSizePolicy(self.Fixed)
            self.gridLayout_advanced.addWidget(self.Picture_advanced_1, 0, 0)

            #advanced right
            self.Picture_advanced_2 = QLabel()
            self.Picture_advanced_2.setAlignment(Qt.AlignCenter)
            self.Picture_advanced_2.setScaledContents(True)
            self.Picture_advanced_2.setSizePolicy(self.Fixed)
            self.gridLayout_advanced.addWidget(self.Picture_advanced_2, 0, 1)






        # Convergence
        self.tab_widget_advanced2= QWidget(self)
        self.gridLayout_advanced2= QGridLayout(self.tab_widget_advanced2)
        self.main_tab.insertTab(i_tab,self.tab_widget_advanced2, 'Convergence')
        self.tables_dict['Convergence'] = i_tab
        self.main_tab.setTabEnabled(i_tab,False)
        i_tab+=1


        #Convergence left
        self.Picture_advanced2_1 = QLabel()
        self.Picture_advanced2_1.setAlignment(Qt.AlignCenter)
        self.Picture_advanced2_1.setScaledContents(True)
        self.Picture_advanced2_1.setSizePolicy(self.Fixed)
        self.gridLayout_advanced2.addWidget(self.Picture_advanced2_1, 0, 0)

        #Convergence right
        self.Picture_advanced2_2 = QLabel()
        self.Picture_advanced2_2.setAlignment(Qt.AlignCenter)
        self.Picture_advanced2_2.setScaledContents(True)
        self.Picture_advanced2_2.setSizePolicy(self.Fixed)
        self.gridLayout_advanced2.addWidget(self.Picture_advanced2_2, 0, 1)



        # Asymmetry        
        
        self.tab_widget_asymmetry= QWidget(self)
        self.verticalLayout_Asymmetry = QVBoxLayout(self.tab_widget_asymmetry)
        self.horizontalLayout_Asymmetry = QHBoxLayout()

        self.main_tab.insertTab(i_tab,self.tab_widget_asymmetry, 'Asymmetry')
        self.main_tab.setTabEnabled(i_tab,False)
        self.tables_dict['Asymmetry'] = i_tab
        i_tab+=1
        
    
        self.Picture_asymmetry = QLabel()
        self.Picture_asymmetry.setAlignment(Qt.AlignCenter)
        self.Picture_asymmetry.setScaledContents(True)
        self.Picture_asymmetry.setSizePolicy(self.Fixed)
        self.horizontalLayout_Asymmetry.addWidget(self.Picture_asymmetry)
        self.verticalLayout_Asymmetry.addLayout(self.horizontalLayout_Asymmetry)

        
        # slider
        self.slider_Asymmetry = QSlider(Qt.Horizontal)
        self.slider_Asymmetry_spin = QSpinBox()
        self.slider_Asymmetry.setVisible(False)
        self.slider_Asymmetry_spin.setVisible(False)
        self.horizontalLayout_Asymmetry2 = QHBoxLayout()
        self.horizontalLayout_Asymmetry2.addWidget(self.slider_Asymmetry)
        self.horizontalLayout_Asymmetry2.addWidget(self.slider_Asymmetry_spin)
        self.verticalLayout_Asymmetry.addLayout(self.horizontalLayout_Asymmetry2)

        # sawtooths         
        self.tab_widget_sawtooths= QWidget(self)
        self.verticalLayout_Sawtooths = QVBoxLayout(self.tab_widget_sawtooths)
        self.horizontalLayout_Sawtooths = QHBoxLayout()

        self.main_tab.insertTab(i_tab,self.tab_widget_sawtooths, 'Sawtooths')
        self.main_tab.setTabEnabled(i_tab,False)
        self.tables_dict['Sawtooths'] = i_tab
        i_tab+=1
        
    
        self.Picture_Sawtooths = QLabel()
        self.Picture_Sawtooths.setAlignment(Qt.AlignCenter)
        self.Picture_Sawtooths.setScaledContents(True)
        self.Picture_Sawtooths.setSizePolicy(self.Fixed)
        self.horizontalLayout_Sawtooths.addWidget(self.Picture_Sawtooths)
        self.verticalLayout_Sawtooths.addLayout(self.horizontalLayout_Sawtooths)


        # Impurities         
        self.tab_widget_impur= QWidget(self)
        self.verticalLayout_Impur = QVBoxLayout(self.tab_widget_impur)
        self.horizontalLayout_Impur = QHBoxLayout()

        self.main_tab.insertTab(i_tab,self.tab_widget_impur, 'Impurities')
        self.main_tab.setTabEnabled(i_tab,False)
        self.tables_dict['Impurities'] = i_tab
        i_tab+=1
        
    
        self.Picture_Impur = QLabel()
        self.Picture_Impur.setAlignment(Qt.AlignCenter)
        self.Picture_Impur.setScaledContents(True)
        self.Picture_Impur.setSizePolicy(self.Fixed)
        self.horizontalLayout_Impur.addWidget(self.Picture_Impur)
        self.verticalLayout_Impur.addLayout(self.horizontalLayout_Impur)

        # slider
        self.slider_Impur = QSlider(Qt.Horizontal)
        self.slider_Impur_spin = QSpinBox()
        self.slider_Impur.setVisible(False)
        self.slider_Impur_spin.setVisible(False)
        self.horizontalLayout_Impur2 = QHBoxLayout()
        self.horizontalLayout_Impur2.addWidget(self.slider_Impur)
        self.horizontalLayout_Impur2.addWidget(self.slider_Impur_spin)
        self.verticalLayout_Impur.addLayout(self.horizontalLayout_Impur2)



        # Poloidal mode number
        self.tab_widget_polnum = QWidget(self)
        self.verticalLayout_polnum = QVBoxLayout(self.tab_widget_polnum)
        self.horizontalLayout_polnum = QHBoxLayout()
        self.main_tab.insertTab(i_tab,self.tab_widget_polnum, 'Poloidal mode num.')
        self.main_tab.setTabEnabled(i_tab,False)
        self.tables_dict['Poloidal mode'] = i_tab
        i_tab+=1

        #wave right
        self.Picture_polnum_1 = QLabel()
        #self.Picture_polnum_1.setText('wave profile left')
        self.Picture_polnum_1.setAlignment(Qt.AlignCenter)
        self.Picture_polnum_1.setScaledContents(True)
        self.Picture_polnum_1.setSizePolicy(self.Fixed)
        self.horizontalLayout_polnum.addWidget(self.Picture_polnum_1)
        #wave left
        self.verticalLayout_polnum.addLayout(self.horizontalLayout_polnum)
        


        i_tab = 0
        # setting tab
        self.left_tab = QTabWidget(self)
        self.gridLayout.addWidget(self.left_tab, 2, 0,2,1)

        # input list
        self.tab_widget_input = QWidget(self)
        self.gridLayout_input = QGridLayout(self.tab_widget_input)
        self.left_tab.insertTab(i_tab,self.tab_widget_input, 'Input')
        self.tables_dict['Input'] = i_tab
        i_tab+=1


        self.Input_label = QLabel("Source:", self.tab_widget_input)
        self.gridLayout_input.addWidget(self.Input_label, 0, 0)
        self.chooseTokamak = QComboBox(self.tab_widget_input)
        item_list = ['Golem - BOLO',
                        "Golem - Camera",
                        "COMPASS - BOLO", 
                        "COMPASS - SXR", 
                        "COMPASS Camera",
                        "JET - slow SXR",
                        "JET - fast SXR" ,
                        "JET - neutrons", 
                        "JET - BOLO", 
                        "ASDEX SXR",
                        "ASDEX SXR_fast",
                        "ASDEX BOLO", 
                        "ASDEX AXUV", 
                        "Tore Supra - SXR",
                        "ToreSupra Camera",
                        "TCV - XTOMO",
                        "TCV - AXUV",
                        "TCV - BOLO",
                        "TCV - DMPX",
                        "TCV - FIR" ,
                        "TCV - XTEPRO",
                        "DIII-D - SXR",
                        "DIII-D - SXR-fast",
                        "DIII-D - BOLO",
                        "DIII-D - DISRAD",
                        ]
        
        
        for item in item_list: 
            self.chooseTokamak.addItem(item)

        self.gridLayout_input.addWidget(self.chooseTokamak, 0, 1)
        self.Shot_label = QLabel("Shot number", self.tab_widget_input)
        self.gridLayout_input.addWidget(self.Shot_label, 1, 0)

        self.lineEdit_Shot = QLineEdit( self.tab_widget_input)
        self.lineEdit_Shot.setMinimumWidth(50)

        self.horizontalLayout_10 = QHBoxLayout()
        self.horizontalLayout_10.addWidget(self.lineEdit_Shot)
        ##add TMIN / TMAX spinboxes
        self.tmin_spin = QDoubleSpinBox()
        self.tmax_spin = QDoubleSpinBox()
        self.tmin_spin.setDecimals(4)
        self.tmax_spin.setDecimals(4)
        self.tmin_spin.setSingleStep(0.1)
        self.tmax_spin.setSingleStep(0.1)


        t_name = self.tokamak.t_name


        self.labelfrom = QLabel( "   Time interval ["+t_name+"]   from:")
        self.labelTo = QLabel( "to")
        self.horizontalLayout_10.addWidget(self.labelfrom)
        self.horizontalLayout_10.addWidget(self.tmin_spin)
        self.horizontalLayout_10.addWidget(self.labelTo)
        self.horizontalLayout_10.addWidget(self.tmax_spin)


        self.gridLayout_input.addLayout(self.horizontalLayout_10, 1, 1)

        self.DataSettingButton = QPushButton("Data Settings", self.tab_widget_input)
        self.gridLayout_input.addWidget(self.DataSettingButton, 2, 1)




        ##output

        self.tab_widget_output = QWidget(self)
        self.gridLayout_output = QGridLayout(self.tab_widget_output)
        self.left_tab.insertTab(i_tab,self.tab_widget_output, 'Output')
        self.tables_dict['Output'] = i_tab
        i_tab+=1

        self.actual_frame = 0

        self.label_output = QLabel("Output Path",self.tab_widget_output)
        self.gridLayout_output.addWidget(self.label_output, 0, 0)

        self.lineEdit_output = QLineEdit(self.tab_widget_output)
        self.gridLayout_output.addWidget(self.lineEdit_output, 0, 1)

        self.AdressButton_output = QToolButton(self.tab_widget_output)
        self.AdressButton_output.setText( "...")
        self.gridLayout_output.addWidget(self.AdressButton_output, 0, 2)


        self.verticalLayout = QVBoxLayout()

        self.enableOutput = QCheckBox("Publication-quality plots", self.tab_widget_output)
        self.PlotAll = QCheckBox("Separated plots", self.tab_widget_output)
        self.Background = QCheckBox("Remove background", self.tab_widget_output)
        self.loglin = QCheckBox("Plot in loglin scale", self.tab_widget_output)

        self.verticalLayout.addWidget(self.enableOutput)
        self.verticalLayout.addWidget(self.PlotAll)
        self.verticalLayout.addWidget(self.Background)
        self.verticalLayout.addWidget(self.loglin)


        self.gridLayout_output.addLayout(self.verticalLayout, 1, 1)

        ##other

        self.tab_widget_other = QWidget(self)
        self.gridLayout_other = QGridLayout(self.tab_widget_other)
        self.left_tab.insertTab(i_tab,self.tab_widget_other, 'Actions')
        self.tables_dict['Actions'] = i_tab
        i_tab+=1

        self.horizontalLayout_act = QHBoxLayout()

        self.verticalLayout_act1 = QVBoxLayout()
        self.verticalLayout_act2 = QVBoxLayout()

        self.reconstruct_box = QCheckBox("Tomography", self.tab_widget_other)
        self.verticalLayout_act1.addWidget(self.reconstruct_box)
        self.SVD_box = QCheckBox("SVD", self.tab_widget_other)
        self.verticalLayout_act1.addWidget(self.SVD_box)
        self.post_box = QCheckBox("Postprocessing", self.tab_widget_other)
        self.verticalLayout_act1.addWidget(self.post_box)
        self.asym_box = QCheckBox("Asymmetries", self.tab_widget_other)
        self.verticalLayout_act1.addWidget(self.asym_box)
        self.saw_box = QCheckBox("Sawtooths", self.tab_widget_other)
        self.verticalLayout_act1.addWidget(self.saw_box)
        self.polnum_box = QCheckBox("Poloidal mode", self.tab_widget_other)
        self.verticalLayout_act2.addWidget(self.polnum_box)
        self.impur_box = QCheckBox("Impurities", self.tab_widget_other)
        self.verticalLayout_act2.addWidget(self.impur_box)
        
        self.gridLayout_other.addLayout(self.verticalLayout_act1, 1, 1)
        self.gridLayout_other.addLayout(self.verticalLayout_act2, 1, 2)

        ## ratio solvers

        self.verticalLayout_solve = QVBoxLayout()

        self.tab_widget_solve = QWidget(self)
        self.gridLayout_solve = QGridLayout(self.tab_widget_solve)
        self.left_tab.insertTab(i_tab,self.tab_widget_solve, 'Solvers')
        self.tables_dict['Solvers'] = i_tab
        i_tab+=1

        self.horizontalLayout_solver_1 = QHBoxLayout()

        self.ratioSolversLabel = QLabel("Ratio presolvers:")
        self.horizontalLayout_solver_1.addWidget(self.ratioSolversLabel)
        self.chooseRatioSolver = QComboBox(self.tab_widget_solve)
        self.chooseRatioSolver.addItem("None")
        self.chooseRatioSolver.addItem("Tikhonov rapid")
        self.chooseRatioSolver.addItem("SVD decomposition")
        self.chooseRatioSolver.addItem("SVD2 decomposition")
        self.chooseRatioSolver.addItem("QR decomposition")
        self.chooseRatioSolver.addItem("GEV decomposition")
        self.chooseRatioSolver.addItem("GSVD decomposition")
        self.chooseRatioSolver.addItem("Tikhonov small resolution")

        #presolves
        self.horizontalLayout_solver_1.addWidget(self.chooseRatioSolver)
        self.verticalLayout_solve.addLayout(self.horizontalLayout_solver_1)

        self.horizontalLayout_solver_2 = QHBoxLayout()

        self.presolverLabel = QLabel("Presolvers:")
        self.horizontalLayout_solver_2.addWidget(self.presolverLabel)
        self.choosePreSolver = QComboBox(self.tab_widget_solve)
        self.choosePreSolver.addItem("None")
        self.choosePreSolver.addItem("Tikhonov rapid")
        self.choosePreSolver.addItem("SVD decomposition")
        self.choosePreSolver.addItem("SVD2 decomposition")
        self.choosePreSolver.addItem("QR decomposition")
        self.choosePreSolver.addItem("GEV decomposition")
        self.choosePreSolver.addItem("GSVD decomposition")
        self.choosePreSolver.addItem("Tikhonov small resolution")

        #solvers

        self.horizontalLayout_solver_2.addWidget(self.choosePreSolver)
        self.verticalLayout_solve.addLayout(self.horizontalLayout_solver_2)

        self.horizontalLayout_solver_3 = QHBoxLayout()
        self.solverLabel = QLabel("Solvers:")
        self.horizontalLayout_solver_3.addWidget(self.solverLabel)
        self.chooseSolver = QComboBox(self.tab_widget_solve)
        self.chooseSolver.addItem("Tikhonov separate")
        self.chooseSolver.addItem("Tikhonov rapid")
        self.chooseSolver.addItem("SVD decomposition")
        self.chooseSolver.addItem("SVD2 decomposition")
        self.chooseSolver.addItem("QR decomposition")
        self.chooseSolver.addItem("GEV decomposition")
        self.chooseSolver.addItem("GSVD decomposition")
        self.chooseSolver.addItem("None")

        self.horizontalLayout_solver_3.addWidget(self.chooseSolver)
        self.verticalLayout_solve.addLayout(self.horizontalLayout_solver_3)
        
        
        #lambda solvers
        self.horizontalLayout_solver_6 = QHBoxLayout()
        self.solverLabel = QLabel("Lambda solvers:")
        self.horizontalLayout_solver_6.addWidget(self.solverLabel)
        self.chooseLamSolver = QComboBox(self.tab_widget_solve)
        self.chooseLamSolver.addItem("GCV")
        self.chooseLamSolver.addItem("PRESS")
        self.chooseLamSolver.addItem("AIC") 
        self.chooseLamSolver.addItem("AICc")    
        self.chooseLamSolver.addItem("BIC")
        self.chooseLamSolver.addItem("CHI2")
        self.chooseLamSolver.addItem("Manual")


        self.horizontalLayout_solver_6.addWidget(self.chooseLamSolver)
        self.verticalLayout_solve.addLayout(self.horizontalLayout_solver_6)
        
        
        
        
        


        self.gridLayout_solve.addLayout(self.verticalLayout_solve, 0, 0)

        self.verticalLayout_solve2 = QVBoxLayout()

        self.spacerItem_s1 = QSpacerItem(40, 20, QSizePolicy.Expanding, QSizePolicy.Minimum)

        if self.tokamak.allow_self_calibration and len(self.tokamak.dets_index) <= 5:
            self.horizontalLayout_solver_4 = QHBoxLayout()
            self.Calb1_label = QLabel("Calibration")
            self.horizontalLayout_solver_4.addItem(self.spacerItem_s1)
            self.horizontalLayout_solver_4.addWidget(self.Calb1_label)
            self.Calb_spins = list()
            for i in range(len(self.tokamak.dets_index)):
                self.Calb_spins.append(QDoubleSpinBox())
                self.Calb_spins[i].setRange(0.01,2)
                self.Calb_spins[i].setSingleStep(0.05)
                self.Calb_spins[i].setDecimals(2)
                self.horizontalLayout_solver_4.addWidget(self.Calb_spins[i])
                

            self.verticalLayout_solve2.addLayout(self.horizontalLayout_solver_4)

        self.verticalLayout_solve2.addItem(self.spacerItem_s1)
        self.horizontalLayout_solver_6 = QHBoxLayout()
        self.Rapid_spin = QSpinBox()
        self.Rapid_label = QLabel("Number of blocks")
        self.Rapid_spin.setRange(1,999)  #one rapid solve is much slower then one non rapid
        self.Rapid_spin.setSingleStep(1)
        
        
        
        self.horizontalLayout_solver_6.addItem(self.spacerItem_s1)
        self.horizontalLayout_solver_6.addWidget(self.Rapid_label)
        self.horizontalLayout_solver_6.addWidget(self.Rapid_spin)
        
        self.horizontalLayout_solver_7 = QHBoxLayout()
        
        self.Fisher_spin = QSpinBox()
        self.Fisher_label = QLabel("Max Fisher steps")
        self.Fisher_spin.setRange(1,20)  #one rapid solve is much slower then one non rapid
       
       
        self.horizontalLayout_solver_7.addItem(self.spacerItem_s1)
        self.horizontalLayout_solver_7.addWidget(self.Fisher_label)
        self.horizontalLayout_solver_7.addWidget(self.Fisher_spin)        
        
        
        self.horizontalLayout_solver_8 = QHBoxLayout()
        self.manual_regularization = QDoubleSpinBox()
        self.manual_reg_lab = QLabel("Regularization")
        self.manual_regularization.setRange(0,1)  
        self.manual_regularization.setSingleStep(0.05)
        self.manual_regularization.setEnabled(False)
  
        self.horizontalLayout_solver_8.addItem(self.spacerItem_s1)
        self.horizontalLayout_solver_8.addWidget(self.manual_reg_lab)
        self.horizontalLayout_solver_8.addWidget(self.manual_regularization)



        
        
        self.horizontalLayout_solver_9 = QHBoxLayout()
        self.horizontalLayout_solver_9.setAlignment(Qt.AlignRight)

        self.Positivity = QCheckBox("Positivity",  self.tab_widget_solve)

        self.horizontalLayout_solver_9.addWidget(self.Positivity)
        
        self.mfi_positivity = QDoubleSpinBox()
        self.mfi_positivity_lab = QLabel(" ")
        self.mfi_positivity.setRange(-1,15)  
        self.mfi_positivity.setSingleStep(0.2)
        self.mfi_positivity.setDecimals(1)
  
        self.horizontalLayout_solver_9.addWidget(self.mfi_positivity_lab)
        self.horizontalLayout_solver_9.addWidget(self.mfi_positivity)



        
        self.verticalLayout_solve2.addLayout(self.horizontalLayout_solver_6)
        self.verticalLayout_solve2.addLayout(self.horizontalLayout_solver_7)
        self.verticalLayout_solve2.addLayout(self.horizontalLayout_solver_8)
        self.verticalLayout_solve2.addLayout(self.horizontalLayout_solver_9)

        

        self.gridLayout_solve.addLayout(self.verticalLayout_solve2, 0,1)



        ################### Transformations !! #############################################
        
        self.verticalLayout_transform = QVBoxLayout()

        self.tab_widget_transform = QWidget(self)
        self.gridLayout_transform = QGridLayout(self.tab_widget_transform)
        self.left_tab.insertTab(i_tab,self.tab_widget_transform, 'Transformations')
        self.tables_dict['Transformations'] = i_tab
        i_tab+=1

        
        self.horizontalLayout_solver_5 = QHBoxLayout()

        self.transform_label = QLabel("Orthogonal transformation ")
        self.chooseTransform = QComboBox(self.tab_widget_transform)
        
        self.transformations = ["None","Abel","Cormack","Fourier-Bessel"]#,"Zoom","Artificial","in-out asymm.",'polynom']
        for t in self.transformations:
            self.chooseTransform.addItem(t)

        self.verticalLayout_transform.addWidget(self.transform_label)
        self.verticalLayout_transform.addWidget(self.chooseTransform)



        horizontalLayout = QHBoxLayout()
        
        spacerItem = QSpacerItem(40, 20, QSizePolicy.Expanding, QSizePolicy.Minimum)
        horizontalLayout.addItem(spacerItem)

        self.transform_label_r = QLabel("Radial order:")
        self.transform_spin_r = QSpinBox()
        self.transform_spin_r.setRange(1,100)

        horizontalLayout.addWidget(self.transform_label_r)
        horizontalLayout.addWidget(self.transform_spin_r)
        
        self.transform_label_a = QLabel("Angular order:")
        self.transform_spin_a = QSpinBox()
        self.transform_spin_a.setRange(0,10)

        horizontalLayout.addWidget(self.transform_label_a)
        horizontalLayout.addWidget(self.transform_spin_a)
        

        self.verticalLayout_transform.addLayout(horizontalLayout)


        

        horizontalLayout = QHBoxLayout()

        spacerItem = QSpacerItem(40, 20, QSizePolicy.Expanding, QSizePolicy.Minimum)
        horizontalLayout.addItem(spacerItem)

        self.divertor_label2 = QLabel("x")


        self.verticalLayout_transform.addLayout(horizontalLayout)


        horizontalLayout = QHBoxLayout()

        spacerItem = QSpacerItem(40, 20, QSizePolicy.Expanding, QSizePolicy.Minimum)


       
        self.gridLayout_transform.addLayout(self.verticalLayout_transform, 0, 0, 1,1)

        self.verticalLayout_transform2 = QVBoxLayout()


        horizontalLayout = QHBoxLayout()

        self.verticalLayout_transform2.addLayout(horizontalLayout)

        
        horizontalLayout = QHBoxLayout()
        horizontalLayout.addItem(spacerItem)

        
        self.boundary_width_label = QLabel("Boundary region width: ")
        self.boundary_width_spin = QDoubleSpinBox()
        self.boundary_width_spin.setRange(0,99)
        self.boundary_width_spin.setDecimals(1)
        self.boundary_width_spin.setSingleStep(0.1)
        horizontalLayout.addWidget(self.boundary_width_label)
        horizontalLayout.addWidget(self.boundary_width_spin)

        self.verticalLayout_transform2.addLayout(horizontalLayout)
         
        
        self.gridLayout_transform.addLayout(self.verticalLayout_transform2,0, 1,1,1)

        self.groupBox_4 = QGroupBox()
        self.gridLayout_4 = QGridLayout(self.groupBox_4)

        self.Allowboundary = QCheckBox("Allow boundary", self.groupBox_4)
        self.gridLayout_4.addWidget(self.Allowboundary, 0, 0)

        self.horizontalLayout_3 = QHBoxLayout()
        self.spacerItem4 = QSpacerItem(40, 20, QSizePolicy.Expanding, QSizePolicy.Minimum)
        self.Precision_label = QLabel("Boundary regularization",self.groupBox_4)
        self.Precision_spin = QDoubleSpinBox(self.groupBox_4)
        self.Precision_spin.setRange(0.0,20.)
        self.Precision_spin.setSingleStep(0.1)
        self.Precision_spin.setDecimals(1)

        self.horizontalLayout_3.addItem(self.spacerItem4)
        self.horizontalLayout_3.addWidget( self.Precision_label )
        self.horizontalLayout_3.addWidget( self.Precision_spin )
        self.gridLayout_4.addLayout(self.horizontalLayout_3, 0, 1)


        self.Smoothing_label = QLabel("Smoothing",self.groupBox_4)
        self.gridLayout_4.addWidget(self.Smoothing_label, 2, 0)

        self.horizontalLayout_4 = QHBoxLayout()
        self.chooseRegularization = QComboBox(self.groupBox_4)
        regularizations = ["Isotropic MFI","Anisotropic MFI","Isotropic MDIFF","Anisotropic MDIFF",
                           "Max Entrophy",'Isotropic MNGR',"Anisotropic MNGR","Anisotropic IMGES"]
        for r in regularizations:
            self.chooseRegularization.addItem(r)
        

        self.spacerItem5 = QSpacerItem(40, 20, QSizePolicy.Expanding, QSizePolicy.Minimum)
        self.horizontalLayout_6 = QHBoxLayout()
        self.Anis_label = QLabel("Anisotropic ratio",self.groupBox_4)
        self.Anis_spin = QDoubleSpinBox(self.groupBox_4)
        self.Anis_spin.setRange(0,99)
        self.Anis_spin.setSingleStep(0.1)
        self.Anis_spin.setDecimals(1)
        self.horizontalLayout_4.addWidget(self.chooseRegularization)
        self.horizontalLayout_4.addItem(self.spacerItem5)
        self.horizontalLayout_4.addWidget(self.Anis_label)
        self.horizontalLayout_4.addWidget(self.Anis_spin)
        self.gridLayout_4.addLayout(self.horizontalLayout_4, 2, 1)

        self.Resolution_label = QLabel("Resolution",self.groupBox_4)
        self.gridLayout_4.addWidget(self.Resolution_label, 3, 0)

        self.nx_spin = QSpinBox( self.groupBox_4)
        self.nx_spin.setRange(10,200)
        self.x_label = QLabel("x",self.groupBox_4)
        self.ny_spin = QSpinBox(self.groupBox_4)
        self.ny_spin.setRange(10,300)
        self.spacerItem3 = QSpacerItem(10, 10, QSizePolicy.Expanding, QSizePolicy.Minimum)
        self.Scale_label = QLabel("Scale errorbars",self.groupBox_4)
        self.Scale_spin = QDoubleSpinBox(self.groupBox_4)
        self.Scale_spin.setRange(0.01,100.0)
        self.Scale_spin.setSingleStep(0.1)
        self.Scale_spin.setDecimals(2)

        self.horizontalLayout_5 = QHBoxLayout()
        self.horizontalLayout_5.addWidget(self.nx_spin)
        self.horizontalLayout_5.addWidget(self.x_label)
        self.horizontalLayout_5.addWidget(self.ny_spin)
        self.horizontalLayout_5.addItem(self.spacerItem3)
        self.horizontalLayout_5.addWidget( self.Scale_label)
        self.horizontalLayout_5.addWidget( self.Scale_spin )

        self.gridLayout_4.addLayout(self.horizontalLayout_5, 3, 1)

        self.gridLayout.addWidget(self.groupBox_4, 2, 1)



        #progress
        self.groupBox_2 = QGroupBox()

        self.gridLayout_5 = QGridLayout(self.groupBox_2)

        self.START_button = QPushButton("START",self.groupBox_2)
        self.gridLayout_5.addWidget(self.START_button, 0, 2)

        self.progressBar = QProgressBar(self.groupBox_2)
        self.gridLayout_5.addWidget(self.progressBar, 0, 0)

        self.gridLayout.addWidget(self.groupBox_2, 3, 1)


        if hasattr(self, 'Calb_spins'):
            for i in range(len(self.tokamak.dets_index)):
                self.Calb_spins[i].setToolTip('Set ratio of detectors.\n!! Can be time dependent!!')
        self.chooseSolver.setToolTip('MFI - solve each timeslice separatly\nMFI rapid - all timeslices together\
            \nfaster but less accurate\nsemianalytical methods \n SVD/QR/GEV/GSVD solver ')
        self.choosePreSolver.setToolTip('Is used to speed up MFI reconstruction')
        self.chooseRatioSolver.setToolTip('Find a relative calibration of the cameras minimizing residuum')
        self.chooseLamSolver.setToolTip('Method to find optimal regularization parameter')
        self.chooseTransform.setToolTip('Select orthogonal transformation applied on the solved matrices')
        self.Allowboundary.setToolTip('The reconstruction will be zero on the boundary')
        self.chooseRegularization.setToolTip('Isotropic (rectangular) smoothing prefers smooth in X and Y direction\nAnisotropic smoothing prefers smooth in magnetic field direction\nMDIFF versions use second derivation instead of first\nif no magnetic data availible, falling back to the isotropic version\nFor anisotropic smoothing is better to use higher resolution >=50')
        self.PlotAll.setToolTip('Export every timeslice in separated plot\nIt can take very long time')
        self.enableOutput.setToolTip('Make a high quality pdf plots, it is slow!!!')

        self.Background.setToolTip('Try to use the first time slices (first 10%) to remove background')
        self.loglin.setToolTip('Transfor data by arcsinh(x) before plotting')

        self.SVD_box.setToolTip('Substract first mode of SVD')
        self.reconstruct_box.setToolTip('Perform tomographic reconstruction')
        self.post_box.setToolTip('Determine center of mass, emissivity,...')
        self.asym_box.setToolTip('Estimate asymmetry of the SXR profile,...')
        self.saw_box.setToolTip('Analyze sawtooth crashes')
        self.impur_box.setToolTip('Estimate impurity density from the radiation profile,...')

        self.polnum_box.setToolTip('Estimate a poloidal number of the mode from the reconstruction')
        self.lineEdit_Shot.setToolTip('Number of shot, expected format is integer, multiple shots are not supported')
        self.Anis_spin.setToolTip('Ratio of anisotropic matrices, B = sigmoid(-x)*Bperp+sigmoid(x)*Bpar \n=> The smaller number the more will reconstruction follow the field\nFor MFI is recommended 0.2, for MDIFF < 0.1' )
        self.Precision_spin.setToolTip('Set the pressure forcing reconstruction to be zero on boundary.\nPrecision of the virtual senzor (boundary) is min(sigma)/X.\nRecommended value is 5 but if it is too unsmooth use lower')
        self.Rapid_spin.setToolTip('Number of equidistant blocks solved using rapid version.\nBe careful, it can create artefacts on the boundaries of the blocks !!')
        self.Fisher_spin.setToolTip('Minimal number of the steps in the Minimum Fisher regularization algorithm')
        self.manual_regularization.setToolTip('Set regularization by hand 0 - smallest, 1 - largest')
        self.Positivity.setToolTip('Force positivity of the solution')
        self.mfi_positivity.setToolTip('Use nonlinear MFI iterations to force positivity')

        self.Scale_spin.setToolTip('Scale the expected errors,\nthe higher number the smoother reconstruction.\nOriginal value is 1')
        self.chooseTokamak.setToolTip('Choose the data source:\nGolem - shots are loaded from internet (shots: 2687-2726)\nCOMPASS - only simulation (shot 2721)\nTore Supra - data from file (shot 2721)\nJET - orig - data from file (shot 65942)\nJET neutrons - data from file (shot 65942)\nJET new - data from file (shots: 65670, 65944 , 65948, 65949, 65952)\n JET bolo - 82291')


        self.slider_rec.valueChanged.connect(self.changeImage_rec)
        self.slider_rec.valueChanged.connect(self.slider_rec_spin.setValue)
        self.slider_rec_spin.valueChanged.connect(self.changeImage_rec)
        self.slider_rec_spin.valueChanged.connect(self.slider_rec.setValue)
        self.slider_SVD.valueChanged.connect(self.changeImage_SVD)
        self.slider_SVD.valueChanged.connect(self.slider_SVD_spin.setValue)
        self.slider_SVD_spin.valueChanged.connect(self.changeImage_SVD)
        self.slider_SVD_spin.valueChanged.connect( self.slider_SVD.setValue)
        self.slider_Asymmetry.valueChanged.connect( self.changeImage_Asym)
        self.slider_Asymmetry.valueChanged.connect(  self.slider_Asymmetry_spin.setValue)
        self.slider_Asymmetry_spin.valueChanged.connect(  self.changeImage_Asym)
        self.slider_Asymmetry_spin.valueChanged.connect(   self.slider_Asymmetry.setValue)


        self.slider_Impur.valueChanged.connect(self.changeImage_Impur)
        self.slider_Impur.valueChanged.connect(self.slider_Impur_spin.setValue)
        self.slider_Impur_spin.valueChanged.connect(self.changeImage_Impur)
        self.slider_Impur_spin.valueChanged.connect(self.slider_Impur.setValue)


        self.AdressButton_output.clicked.connect( self.showDialogOutput)
        self.DataSettingButton.clicked.connect( self.showSetting)
        self.START_button.clicked.connect( self.actionStart)
        self.chooseRegularization.currentIndexChanged.connect(self.regularizationChanged)
        self.Allowboundary.stateChanged.connect(self.Precision_label.setEnabled)
        self.Allowboundary.stateChanged.connect(self.Precision_spin.setEnabled)

        self.enableOutput.stateChanged.connect(self.lineEdit_output.setEnabled)
        self.enableOutput.stateChanged.connect(self.AdressButton_output.setEnabled)
        self.enableOutput.stateChanged.connect(self.enableOutputChanged)

        self.SVD_box.stateChanged.connect(self.toggleSVD)
        self.reconstruct_box.stateChanged.connect(self.toggle_reconstruct)
        self.post_box.stateChanged.connect(self.toggle_post)
        self.asym_box.stateChanged.connect(self.toggle_asym)
        self.saw_box.stateChanged.connect(self.toggle_saw)
        self.polnum_box.stateChanged.connect(self.toggle_polnum)
        self.impur_box.stateChanged.connect(self.toggle_impur)


        self.chooseRatioSolver.currentIndexChanged.connect(self.toggle_calb)
        self.chooseLamSolver.currentIndexChanged.connect(self.allow_LamSolv)
        self.chooseTransform.currentIndexChanged.connect(self.toggle_transform)
        self.chooseSolver.currentIndexChanged.connect(self.toggle_presolver)
        self.Background.stateChanged.connect(self.toggle_presolver)

        self.lineEdit_Shot.editingFinished.connect( self.change_shot)
        self.main_tab.currentChanged.connect(self.updatePanels)




        self.thread = MainThread(self)
        self.thread.output.connect(self.updateUi)
        self.thread.failed.connect(self.failed)
        self.thread.newValue.connect( self.progressBar.setValue)
        self.thread.newValue.connect( self.progressBar.setValue)
        self.thread.finished.connect( self.changeStartButton)


        self.setCentralWidget(self.centralwidget)
        self.resize(self.setting['screen_width']/1.5, self.setting['screen_height']/1.5)
        self.setCenter()

        self.setValues()

        self.tmin_spin.valueChanged.connect( self.tminChanged)
        self.tmax_spin.valueChanged.connect( self.tmaxChanged)
        


    def toggle_calb(self, value):
        value = ( value == 0 )
        if self.tokamak.allow_self_calibration and hasattr(self,'Calb_spins'):
            for i in range(len(self.tokamak.dets_index)):
                    self.Calb_spins[i].setEnabled(value)
    def allow_LamSolv(self, value):
        self.manual_regularization.setEnabled( value == 6) #manual


    def toggle_presolver(self, value):
        value = self.chooseSolver.currentIndex()==0 or self.Background.isChecked()
        self.choosePreSolver.setEnabled(value)
        if self.chooseSolver.currentIndex() in [0,1,7]:
            self.manual_regularization.setEnabled(False)
            self.chooseLamSolver.setEnabled(False)
        else:
            self.chooseLamSolver.setEnabled(True)
            self.allow_LamSolv(self.chooseLamSolver.currentIndex() )


    def toggle_reconstruct(self, value):
        self.SVD_box.setEnabled(value)
        self.post_box.setEnabled(value)
        self.groupBox_4.setEnabled(value)
        self.post_box.setEnabled(value)
        self.asym_box.setEnabled(value)
        self.saw_box.setEnabled(value)
        self.impur_box.setEnabled(value)
        self.polnum_box.setEnabled(value)


        self.main_tab.setTabEnabled(self.tables_dict['Main'],value)
        self.left_tab.setTabEnabled(self.tables_dict['Solvers'],value)
        self.left_tab.setTabEnabled(self.tables_dict['Transformations'],value)



    def toggleSVD(self, value):
        self.main_tab.setTabEnabled(self.tables_dict['SVD'],value)



    def toggle_post(self,value):
        self.main_tab.setTabEnabled(self.tables_dict['Position'],value)
        self.main_tab.setTabEnabled(self.tables_dict['Emissivity'],value)
        if 'Equilibrium' in self.tables_dict:
            self.main_tab.setTabEnabled(self.tables_dict['Equilibrium'],value)
        self.main_tab.setTabEnabled(self.tables_dict['Convergence'],value)
    def toggle_asym(self,value):
        self.main_tab.setTabEnabled(self.tables_dict['Asymmetry'],value)

    def toggle_saw(self,value):
        self.main_tab.setTabEnabled(self.tables_dict['Sawtooths'],value)

    def toggle_impur(self,value):
        self.main_tab.setTabEnabled(self.tables_dict['Impurities'],value)

    def toggle_polnum(self,value):
        self.main_tab.setTabEnabled(self.tables_dict['Poloidal mode'],value)


    def toggle_transform(self,value):

        self.transform_label_r.setEnabled(value in [1,2,3,5,6,7])
        self.transform_spin_r.setEnabled(value in [1,2,3,5,6,7] )
        self.transform_label_a.setEnabled(value in [1,2,3,5,6,7])
        self.transform_spin_a.setEnabled(value in [1,2,3,5,6,7] )

    def toggle_boundary(self, value):
        self.boundary_width_label.setEnabled( not value)
        self.boundary_width_spin.setEnabled( not value)
        self.boundary_width_spin.setEnabled( not value)
        

    def enableOutputChanged(self, value):
        pass


    def setCenter(self):
        frect = QDesktopWidget.frameGeometry(self)
        frect.moveCenter(QDesktopWidget().availableGeometry(self).center());
        self.move(frect.topLeft())

    def showSetting(self):

        self.getValues()

        setting, tokamak =  TestData(self)

        if setting is None:
            return
        self.setting = setting
        self.tokamak = tokamak
        self.tokamak2 = tokamak

        from data_preprocessing import DataSettingWindow
        win=DataSettingWindow(self)
        win.show()
        self.setValues()

    def showDialogOutput(self):
        if os.path.isdir(self.lineEdit_output.text()):
            path = self.lineEdit_output.text()
        else:
            path = os.path.normpath(os.path.dirname(sys.argv[0]))
        DirectoryNameOutput = QFileDialog.getExistingDirectory(self, 'Save to Directory', path)
        if DirectoryNameOutput is not None:
            self.lineEdit_output.setText(str(DirectoryNameOutput))

    def regularizationChanged(self, value):
        if value in (0,2,4,5):
            state = False
        else:
            state=True
        self.Anis_spin.setEnabled(state)
        self.Anis_label.setEnabled(state)
        
 
        
    def tminChanged(self,value):
        if self.tmax < value:
            self.tmax_spin.setValue(value)
        self.tmin_spin.setValue(value)
        self.tmin = value
        self.tokamak.tmin = self.tmin

    def tmaxChanged(self,value):
        if self.tmin > value:
            self.tmin_spin.setValue(value)
        self.tmax_spin.setValue(value)
        self.tmax = value
        self.tokamak.tmax = self.tmax

    def getValues(self):
        #print 'getValues',self.tokamak.index
        """
        Load data from GUI and save them to file
        """
 

        self.tmax = self.tmax_spin.value()
        self.tmin = self.tmin_spin.value()
        self.tokamak.tmax = self.tmax
        self.tokamak.tmin = self.tmin

        self.gui = True
        try:
            self.shot = int(self.lineEdit_Shot.text())
        except:
            QMessageBox.warning(self,"Wrong discharge number", "Write here an appropriate discharge number ",QMessageBox.Ok)

        self.setting['shot'] = self.shot

        self.rapid_blocks = self.Rapid_spin.value()
        self.ifishmax = self.Fisher_spin.value()
        self.man_regularization = self.manual_regularization.value()

        self.solver =self.chooseSolver.currentIndex()
        self.presolver =self.choosePreSolver.currentIndex()
        self.ratiosolver = self.chooseRatioSolver.currentIndex()
        self.lambda_solver = self.chooseLamSolver.currentIndex()
        self.lambda_solver = ['gcv','press','aic','aicc','bic','chi2','manual'][self.lambda_solver]
        if self.lambda_solver == 'manual':
            self.lambda_solver = str(self.man_regularization)
 
        if self.ratiosolver == 0 and self.tokamak.allow_self_calibration:
            config.calb = ones(self.tokamak.Ndets)
            for i in range(self.tokamak.Ndets):
                try:
                    config.calb[i] = self.Calb_spins[i].value()
                except:
                    pass
                
                
        self.nx=self.nx_spin.value()
        self.ny=self.ny_spin.value()
        self.regularization= self.chooseRegularization.currentIndex()
        self.plot_all = self.PlotAll.isChecked()
        self.rem_back = self.Background.isChecked()
        self.plot_loglin = self.loglin.isChecked()

        self.plot_svd = self.SVD_box.isChecked()
        self.reconstruct = self.reconstruct_box.isChecked()
        self.plot_poloidal_spectrum = self.polnum_box.isChecked()
        self.postprocessing = self.post_box.isChecked()
        self.asymmetry = self.asym_box.isChecked()
        self.sawtooths = self.saw_box.isChecked()
        self.impurities = self.impur_box.isChecked()
        self.postprocessing |= self.impurities
        self.tok_index = self.chooseTokamak.currentIndex()
        self.enable_output = self.enableOutput.isChecked()
        self.positive_constrain = self.Positivity.isChecked()
        self.rgmin = 10**(-self.mfi_positivity.value())

      
        if self.enable_output:
            self.output_path = str(self.lineEdit_output.text())

        self.danis= self.Anis_spin.value()
        self.boundary = self.Precision_spin.value() if self.Allowboundary.isChecked() else -1
        self.error_scale = self.Scale_spin.value()


        #==========  transformations =========================
        self.boundary_width = self.boundary_width_spin.value()
        self.transform_order_a = self.transform_spin_a.value()
        self.transform_order_r = self.transform_spin_r.value()

        self.transform_index = self.chooseTransform.currentIndex()


        #BUG very ugly solution!!
        for a in list(vars(self.tokamak).keys()):              # load all data from class
            if a not in ('shot','tok_index') and hasattr(self,a):
                setattr(self.tokamak,a, getattr(self, a))

        for i in self.setting:
            if hasattr(self,i):
                self.setting[i] = getattr(self,i)


    def setValues(self):
        """
        Load data from dictionary and set them to the GUI
        """
        
        for i in self.setting:
            setattr(self,i, self.setting[i])
        
        
        try:
            self.tokamak = self.tokamak_tmp
        except:
            pass

        for a in list(vars(self.tokamak).keys()):              # load all data from class
            setattr(self,a, getattr(self.tokamak,a))
            
        self.Allowboundary.setChecked(self.boundary>=0)
        self.Rapid_spin.setValue(self.rapid_blocks)
        self.PlotAll.setChecked(self.plot_all)
   
        self.Fisher_spin.setValue(self.ifishmax)

        self.Background.setChecked(self.rem_back)
        self.loglin.setChecked(self.plot_loglin)

        self.nx_spin.setValue(int(self.nx))
        self.ny_spin.setValue(int(self.ny))
        self.lineEdit_output.setText(self.output_path)
        self.lineEdit_Shot.setText(str(self.shot))
        self.chooseTokamak.setCurrentIndex(self.tok_index)
        self.chooseRegularization.setCurrentIndex(self.regularization)

        self.chooseSolver.setCurrentIndex(self.solver)
        self.choosePreSolver.setCurrentIndex(self.presolver)
        self.chooseRatioSolver.setCurrentIndex(self.ratiosolver)
        
        lam_solvers = ['gcv','press','aic', 'aicc','bic','chi2','manual']
        if self.lambda_solver in lam_solvers:
            ind = [i for i,n in enumerate( lam_solvers) if n == self.lambda_solver ]
            self.chooseLamSolver.setCurrentIndex(ind[0])
            self.manual_regularization.setValue(0.75)
        else:
            try: self.manual_regularization.setValue(float(self.lambda_solver))
            except: raise Exception( 'unknown lambda solver '+str(self.lambda_solver))
            self.chooseLamSolver.setCurrentIndex(6) #manual

        self.Positivity.setChecked(self.positive_constrain)
        self.mfi_positivity.setValue(-log10(self.rgmin))


        if hasattr(self,'Calb_spins'):
            if all([ isscalar(i) for i in self.tokamak.calb_0]):
                calb =         self.tokamak.get_calb()
                for i in range(self.tokamak.Ndets):
                    self.Calb_spins[i].setValue(calb[i])


        if self.regularization in (0,2,4,6):
            self.Anis_spin.setEnabled(False)
            self.Anis_label.setEnabled(False)
        self.enableOutput.setChecked(self.enable_output)

        if self.enableOutput == False:
            self.lineEdit_output.setEnabled(False)
            self.AdressButton_output.setEnabled(False)
    
        self.SVD_box.setChecked(self.plot_svd)
        self.post_box.setChecked(self.postprocessing)
        self.asym_box.setChecked(self.asymmetry)
        self.saw_box.setChecked(self.sawtooths)
        self.impur_box.setChecked(self.impurities)

        self.polnum_box.setChecked(self.plot_poloidal_spectrum)
        self.reconstruct_box.setChecked(self.reconstruct)

        self.toggle_post(self.postprocessing)
        self.toggle_asym(self.asymmetry)
        self.toggle_saw(self.sawtooths)
        self.toggle_impur(self.impurities)

        self.toggleSVD(self.plot_svd)
        self.toggle_polnum(self.plot_poloidal_spectrum)
        self.toggle_reconstruct(self.reconstruct)


        self.Precision_spin.setEnabled(self.boundary>=0)
        self.Precision_label.setEnabled(self.boundary>=0)
        self.Anis_spin.setValue(self.danis)
        self.Precision_spin.setValue(double(self.boundary))
        self.Scale_spin.setValue(self.error_scale)


        #==========  transformations =========================
        self.boundary_width_spin.setValue(self.boundary_width)
        self.transform_spin_r.setValue(self.transform_order_r)
        self.transform_spin_a.setValue(self.transform_order_a)

        self.chooseTransform.setCurrentIndex(-1)
        self.chooseTransform.setCurrentIndex(self.transform_index)
        
        self.toggle_transform(self.transform_index)

        self.tmin_spin.setRange(double(self.tokamak.min_tvec), double(self.tokamak.max_tvec))
        self.tmax_spin.setRange(double(self.tokamak.min_tvec), double(self.tokamak.max_tvec))

        self.tmin_spin.setValue(double(self.tmin))
        self.tmax_spin.setValue(double(self.tmax))

    def changeStartButton(self, state=0):
        if state == 0:
            self.START_button.setText('START')
            try:
                self.START_button.clicked.disconnect( self.killThread)
            except:
                pass
            self.START_button.clicked.connect( self.actionStart)
        else:
            self.START_button.setText('STOP')
            self.START_button.clicked.disconnect( self.actionStart)
            self.START_button.clicked.connect(self.killThread)


    def killThread(self):
        """
        Stop thread when the STOP button as pressed. Sometimes it must waits to end of last cycle
        """
        if self.thread.isRunning():
            print("\n**********KILL********** (waits to end of one cycle)\n")
            self.thread.requestInterruption()
            self.thread.wait()
   
        self.thread.output.disconnect( self.updateUi)
        self.thread.failed.disconnect( self.failed)
        self.thread.newValue.disconnect( self.progressBar.setValue)
        self.thread.finished.disconnect( self.changeStartButton)
        try:
            self.START_button.clicked.disconnect( self.killThread)
        except:
            self.START_button.clicked.disconnect( self.actionStart)
        del self.thread                #there were problems when the thread is killed during multithread solve

        self.thread = MainThread(self)

        self.changeStartButton(0)

        self.thread.output.connect(self.updateUi)
        self.thread.failed.connect(self.failed)
        self.thread.newValue.connect(self.progressBar.setValue)
        self.thread.finished.connect(self.changeStartButton)
        self.progressBar.reset()


    def updateUi(self, output_list = []):
        """
        Refresh previews in GUI
        """

        self.progress = self.thread.progress   # used to control progressBar !!!


        self.actual_frame = 0


        if self.reconstruct:
            tsteps = self.setting.get('tsteps',2)

        self.changeImage_rec(self.actual_frame)
        if self.plot_all and self.reconstruct and self.tmin!= self.tmax:
            self.slider_rec.setVisible(True)
            self.slider_rec.setValue(self.actual_frame)
            self.slider_rec_spin.setVisible(True)
            self.slider_rec_spin.setMaximum(tsteps)
            self.slider_rec.setMaximum(tsteps)
        else:
            self.slider_rec.setVisible(False)
            self.slider_rec_spin.setVisible(False)





        if self.plot_svd and self.reconstruct and len(output_list):
            try:
                make_svd(*output_list)
            except Exception as e:
                QMessageBox.warning(self,"SVD problem", "SVD  failed\n "+str(e),QMessageBox.Ok)
                return 1


            self.changeImage_SVD(self.actual_frame)


        if self.plot_all and self.reconstruct and tsteps > 1:
            self.slider_SVD.setVisible(True)
            self.slider_SVD.setValue(self.actual_frame)
            self.slider_SVD_spin.setVisible(True)
            self.slider_SVD_spin.setMaximum(tsteps)
            self.slider_SVD.setMaximum(tsteps)
        else:
            self.slider_SVD.setVisible(False)
            self.slider_SVD_spin.setVisible(False)
            
        if self.plot_all and self.reconstruct and tsteps > 1:
            self.slider_Asymmetry.setVisible(True)
            self.slider_Asymmetry.setValue(self.actual_frame)
            self.slider_Asymmetry_spin.setVisible(True)
            self.slider_Asymmetry_spin.setMaximum(tsteps)
            self.slider_Asymmetry.setMaximum(tsteps)
            self.changeImage_Asym(self.actual_frame)

        else:
            self.slider_Asymmetry.setVisible(False)
            self.slider_Asymmetry_spin.setVisible(False)
            self.changeImage_Asym(self.actual_frame)
                
            
        
        if self.postprocessing:
            images = {  self.Picture_position_1:'xmass',\
                        self.Picture_position_2:'ymass',\
                        self.Picture_emis_1:'power',\
                        self.Picture_emis_2:'profile',\

                        self.Picture_advanced2_1:'lam',\
                        self.Picture_advanced2_2:'chi2'  }
            if self.post_proces_equi:
                images.update({self.Picture_advanced_1:'shafr_shift',\
                        self.Picture_advanced_2:'elongation'})
            
            for image, name in list(images.items()):
                img = QPixmap(self.tmp_folder+'%s_%d.png'%(name,self.shot))
                img.setDevicePixelRatio(qtscale)
                image.setPixmap(img)
                
        if self.sawtooths:
            image = QPixmap(self.tmp_folder+'sawtooths%d.png'%self.shot)
            image.setDevicePixelRatio(qtscale)
            self.Picture_Sawtooths.setPixmap(image)

        if self.plot_all and self.reconstruct and tsteps > 1:
            self.slider_Impur.setVisible(True)
            self.slider_Impur.setValue(self.actual_frame)
            self.slider_Impur_spin.setVisible(True)
            self.slider_Impur_spin.setMaximum(tsteps)
            self.slider_Impur.setMaximum(tsteps)
            self.changeImage_Impur(self.actual_frame)

        else:
            self.slider_Impur.setVisible(False)
            self.slider_Impur_spin.setVisible(False)
            self.changeImage_Impur(self.actual_frame)
            

        if self.plot_poloidal_spectrum:
            prewA = self.tmp_folder+'PoloidalModeNumber.png'
            image = QPixmap(prewA)
            image.setDevicePixelRatio(qtscale)
            self.Picture_polnum_1.setPixmap(image)
            self.main_tab.setCurrentIndex(5)
     
        self.progressBar.setValue(100)
        self.killThread()


        print('DONE')

    def changeImage_rec(self, number):
        self.actual_frame = number

        if number == 0:
            prewA = self.tmp_folder+'/previewA%d_rec_.png'%self.shot
            prewB = self.tmp_folder+'/previewB%d_rec_.png'%self.shot
        else:
            name = '%d_%.4d'%(self.shot,number-1)
            prewA = self.tmp_folder+'/brightness_'+name+'_rec_.png'
            prewB = self.tmp_folder+'/emissivity_'+name+'_rec_.png'


        pix_map = QPixmap(prewA) 
        pix_map.setDevicePixelRatio(qtscale)
        self.Picture_recons_1.setPixmap(pix_map)
        pix_map = QPixmap(prewB)
        pix_map.setDevicePixelRatio(qtscale)
        self.Picture_recons_2.setPixmap(pix_map)

    def change_shot(self):
        self.getValues()
        self.change_overview_plots()

    def change_overview_plots(self):
    
        img_path = self.local_path+'/geometry/'+self.tokamak.name+'/overview_plots/overview_%d.png'%self.setting['shot']
        if os.path.isfile(img_path):
            image = QPixmap(img_path)
            image.setDevicePixelRatio(qtscale)
            rect = QRect(78, 60, 730, 485)
            image = image.copy(rect)
            self.OverviewPlot.setPixmap(image)

    def changeImage_SVD(self, number):
        self.actual_frame = number

        if number == 0:
            prewA = self.tmp_folder+'/SVD_%d.png'%self.shot
            prewB = '' 
        else:
            name = '%d_%.4d'%(self.shot,number-1)
            prewA = self.tmp_folder+'/brightness_'+name+'_SVD_.png'
            prewB = self.tmp_folder+'/emissivity_'+name+'_SVD_.png'
        
        image = QPixmap(prewA)
        image.setDevicePixelRatio(qtscale)
        self.Picture_SVD_1.setPixmap(image)
        
        image = QPixmap(prewB)
        image.setDevicePixelRatio(qtscale)
        self.Picture_SVD_2.setPixmap(prewB)        
        
    def updatePanels(self,panel_ind):
        #update position in the actual panel
        name = [name for name,ind in list(self.tables_dict.items()) if ind == panel_ind]

        if 'Main' in name: 
            self.changeImage_rec(self.actual_frame)
            self.slider_rec.setValue(self.actual_frame)
        if 'SVD' in name: 
            self.changeImage_SVD(self.actual_frame)
            self.slider_SVD.setValue(self.actual_frame)
        if 'Asymmetry'in name: 
            self.changeImage_Asym(self.actual_frame)
            self.slider_Asymmetry.setValue(self.actual_frame)
        if 'Impurities' in name:  
            self.changeImage_Impur(self.actual_frame)
            self.slider_Impur.setValue(self.actual_frame)

        
    def changeImage_Asym(self, number):
        self.actual_frame = number

        if number == 0 and 'tsteps' in self.setting and self.setting['tsteps'] > 1:
            prew =  self.tmp_folder+'/asymmetry_2D_%d.png'%self.shot
        elif 'tsteps' in self.setting and self.setting['tsteps'] > 1:
            prew =  self.tmp_folder+'/asymmetry_plot_%.4d.png'%(number-1)
        else:
            prew =  self.tmp_folder+'/asymmetry_plot_%.4d.png'%0
            
            
        image = QPixmap(prew)
        image.setDevicePixelRatio(qtscale)
        self.Picture_asymmetry.setPixmap(image)
        
        
    def changeImage_Impur(self, number):
        self.actual_frame = number

        if number == 0 and 'tsteps' in self.setting and self.setting['tsteps'] > 1:
            prew =  self.tmp_folder+'/impurities_2D_%d.png'%self.shot
        elif 'tsteps' in self.setting and self.setting['tsteps'] > 1:
            prew =  self.tmp_folder+'/imp_plot_%.4d.png'%(number-1)
        else:
            prew =  self.tmp_folder+'/imp_plot_%.4d.png'%0
        
        image = QPixmap(prew)
        image.setDevicePixelRatio(qtscale)
        self.Picture_Impur.setPixmap(image)
        

      
    


    def failed(self, error):
        print("Reconstruction failed ")
        gc.collect()

        if error == 2:
            QMessageBox.warning(self ,"Error",  "Not enough of operation memory,\
                            select shorter time interval" ,QMessageBox.Ok)
        else:
            QMessageBox.warning(self,"Execution problem", "Reconstruction failed ",QMessageBox.Ok)

    def actionStart(self):
        """
        Function called when START button pressed.
        """
        print('START...')

        self.getValues()

        setting, tokamak =  TestData(self)
        if setting is None:
            return
        

        self.setting = setting
        self.tokamak = tokamak


        ##########    setup tokamak ######################

        self.tokamak.prepare_tokamak()
        self.progressBar.setValue(1)

        ind_min,ind_max = self.tokamak.tvec.searchsorted([self.tmin,self.tmax])+1
        steps = int((self.tmax-self.tmin)*self.tokamak.sample_freq)//self.data_undersampling
        
     
        ncpu = 1
        try:
            from multiprocessing import cpu_count
            ncpu = min(8,cpu_count())
        except:
            pass

        if rcParams['backend'].lower() == 'agg':
            speed = 100*ncpu
        else:
            speed = 100
        if self.enableOutput.isChecked():
            speed /= 10
            
        
        if (steps > speed) and self.plot_all and self.reconstruct:
            reply = QMessageBox.question(self, 'Message', "Are you sure that you want to plot all "+str(steps)
                                        +" snapshots? \n It can take quite long...", QMessageBox.Yes, QMessageBox.No)
            if reply == QMessageBox.No:
                self.progressBar.reset()
                return
        
        try:
            FreeMemory=int(os.popen("free -m").readlines()[2].split()[3])
        except:
            FreeMemory = 1e4
        if (double(steps)*double(self.nx*self.ny)*16.0  > FreeMemory*1e6 and self.reconstruct) and  not (self.solver == 2 and self.ratiosolver in [0,2]):   #limit of memory usage
            print('free memory ', FreeMemory*1e6, double(steps)*double(self.nx*self.ny)*16.0)
            QMessageBox.warning(self,"Setting error", "Array is too big \n Use lower resolution or less timeslices"  ,QMessageBox.Ok)
            return
      
      
        import shutil


        try:
            
            debug('clean tmp')
            #print('BUG no clean tmp')
            shutil.rmtree(setting['tmp_folder'])
        except Exception as e:
            print('cleaning of the TMP folder %s has failured'%setting['tmp_folder'], e)
            print(os.listdir(setting['tmp_folder'] ))
        
        
        try:
            if not os.path.exists(setting['tmp_folder']):
                os.mkdir(setting['tmp_folder'])
        except Exception as e:
            time.sleep(0.1)
            os.mkdir(setting['tmp_folder'])
                
 
       
        if self.reconstruct:
            self.thread.prepare(self.setting, self.tokamak)
            self.thread.start()
            self.changeStartButton(1)
        elif self.plot_poloidal_spectrum  and not self.reconstruct :
            self.updateUi(None)
        else:
            QMessageBox.warning(self,"Missing action", 'Select an action on "Actions" list' ,QMessageBox.Ok)

        try:
            from sksparse.cholmod import  cholesky
            
        except:
            QMessageBox.warning(self,"Import Error", "Can't find Cholesky decomposition for sparse matrices (module scikit.sparse) !!! \n Falling back to slower version" ,QMessageBox.Ok)

        self.resize(self.setting['screen_width']/1.5, self.setting['screen_height']/1.5)


class MainThread(QThread):
    """
    Class running  in separated thread, call main solving function "main.py"
    """
    newValue = pyqtSignal(int)
    failed = pyqtSignal(int)
    output = pyqtSignal(tuple)

    
    
    def __init__(self, parent):
        QThread.__init__(self, parent)
        self.progress = IterateProgress(self,parent.progressBar)
        self.parent = parent
    def prepare(self,setting, tokamak):
        self.setting = setting
        self.tokamak = tokamak
    def run(self):
        print('RUN')

        try:
            output_list = tomography(self.setting, self.tokamak,  self.progress)
        except KeyboardInterrupt:
            self.newValue.emit(0)
            return
        
        except MemoryError as e:
            print("memory error detected"+ str(e))
            self.failed.emit(2)
            self.newValue.emit(0)
            
            self.failed.emit( 2)
            self.newValue.emit(0)
            
            
            return
        except Exception as e:
            print("some error detected"+ str(e))
            self.failed.emit(1)
            self.newValue.emit(0)
            traceback.print_exc()
            return
        
        if hasattr(self,'isInterruptionRequested') and self.isInterruptionRequested():
            self.newValue.emit(0)
            return
        
        try:
            make_graphs(output_list, False)        # gnuplot can plot in separated thread
            self.progress.iterateStep()
            if self.parent.plot_svd:
                make_graphs(copy(output_list), True)
                self.progress.iterateStep()


        except KeyboardInterrupt:
            self.newValue.emit(0)
        except Exception as e:
            print("Plotting failed"+ str(e))
            self.failed.emit( 1)
            self.newValue.emit(0)
            traceback.print_exc()
            return
        
        
        
        if hasattr(self,'isInterruptionRequested') and self.isInterruptionRequested():
            self.newValue.emit(0)
            traceback.print_exc()
            return
        
        try:

            
            if self.setting['postprocessing'] and self.setting['reconstruct'] and self.setting['tsteps']>1:
                postprocessing_plot(output_list)
                self.progress.iterateStep()

            if self.setting['sawtooths'] and self.setting['reconstruct']:
                sawtooths_detection(output_list)
                self.progress.iterateStep()

            if self.setting['impurities'] and self.setting['reconstruct']:
                imp_analysis(output_list)
                self.progress.iterateStep()

            if self.setting['plot_poloidal_spectrum'] and self.setting['reconstruct']:
                CalcPoloidalModeSpectrum(output_list)
                self.progress.iterateStep()
            
            if self.setting['asymmetry'] and self.setting['reconstruct']:
                from asymmetries import  EvalAsymmetry
                EvalAsymmetry(  output_list)
                self.progress.iterateStep()

        except KeyboardInterrupt:
            self.newValue.emit(0)
        except Exception as e:
            print("Asymmetry failed"+ str(e))
            self.failed.emit( 1)
            self.newValue.emit(0)
            traceback.print_exc()
            return
   
        self.output.emit(output_list)


class IterateProgress(QObject):
    """
    Class controling the progress bar
    """
     
    
    def __init__(self,parent,progressBar):
        self.progressBar = progressBar
        self.parent = parent
        self.value = 0
        self.Step = 0
        
    def setNumSteps(self, num):
        self.numSteps = num
        
    def setNumSubSteps(self, num):
        self.numSubSteps = num
        
    def iterateStep(self, steps = 1):
        self.Step += steps
        self.value = max([self.Step * 100.0/(self.numSteps), self.value])
        self.parent.newValue.emit(int(self.value))
        
    def iterateSubStep(self):
        self.value += 100.0/(self.numSteps*self.numSubSteps)
        self.parent.newValue.emit(int(self.value))
        
    def setValue(self,value):
        self.value = value
        self.parent.newValue.emit(int(self.value))


def TestData(tomography_object):
    """
    Check if data settings are correct
    """

    tomography_object.getValues()
    

    new_setting = 1
    setting = tomography_object.setting
    tokamak = tomography_object.tokamak
    setting['tokamak_tmp'] = tokamak
    
    
    

    try:
        setting['shot'] = int(tomography_object.lineEdit_Shot.text())  #shot is expected as integer, multiple shots are not supported
        setting['nx']=tomography_object.nx_spin.value()
        setting['ny']=tomography_object.ny_spin.value()

        if setting['tok_index'] != tomography_object.chooseTokamak.currentIndex():
            tokamak.wrong_dets = []


        setting['tmin']=tomography_object.tmin
        setting['tmax']=tomography_object.tmax

        enable_output = tomography_object.enableOutput.isChecked()
        if enable_output:
            try:
                 output_path = str(tomography_object.lineEdit_output.text())
            except:
                QMessageBox.warning(tomography_object,"Input error", "Output path couldn't be loaded \n Disabling output"  ,QMessageBox.Ok)
                enable_output = False
                tomography_object.enableOutput.setChecked(False)
                tomography_object.lineEdit_output.setEnabled(False)
                tomography_object.AdressButton_output.setEnabled(False)
                return

    except NameError as details:
        QMessageBox.warning(tomography_object,"Input error", details.message ,QMessageBox.Ok)
        return None,None
    except:
        QMessageBox.warning(tomography_object,"Input error", "Check input values, shot is expected as integer"  ,QMessageBox.Ok)
        return None,None
    if enable_output:
        if not os.path.isdir( output_path):
            os.makedirs(output_path)

    try:
        new_setting = setting.copy()
        new_tokamak = loaddata(setting)
        new_setting['tokamak_tmp'] = new_tokamak
        

        if enable_output:
            os.listdir(output_path)
    except IOError as details:
        QMessageBox.warning(tomography_object,"Error", "File doesn't exist: "+ str(details.filename)   ,QMessageBox.Ok)
        #raise
        return None,None
    except OSError as details:
        QMessageBox.warning(tomography_object,"Error", "Folder doesn't exist: "+ str(details.filename)   ,QMessageBox.Ok)
        #raise
        return None,None
    except NameError as details:
        QMessageBox.warning(tomography_object,"Error", details.message ,QMessageBox.Ok)
        #raise
        return None,None
    except MemoryError:
        QMessageBox.warning(tomography_object,"Error",  "Not enough of operation memory, select shorter time interval" ,QMessageBox.Ok)
        #raise
        return None,None
    except:
        QMessageBox.warning(tomography_object,"Unknown Error", "Data couldn't be loaded\nIf you changed tokamak now\ntry to remove  tomography.npy file",QMessageBox.Ok)
        #raise
        return None,None
    List =  os.listdir('./')

    return new_setting, new_tokamak   #!!! a teď není potřeba aby se to načítalo 2x !!!!!!!!!!

def main(setting, tokamak):
    


    app = QApplication(sys.argv)
    app.setStyle(QStyleFactory.create("plastique"))
    new_font = app.font()
    new_font.setPointSize(  12 )
    app.setFont( new_font )
    setting['screen_width']  = app.desktop().screenGeometry().width()
    setting['screen_height'] = app.desktop().screenGeometry().height()

    myapp = MainWindow(setting, tokamak)
    myapp.show()
    app.exec_()



