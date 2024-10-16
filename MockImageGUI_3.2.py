import yt
import os
import sys
import cmyt
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PyQt5.QtWidgets import QApplication, QWidget, QLabel, QLineEdit, QPushButton, QVBoxLayout, QHBoxLayout, QFileDialog, QCheckBox, QStackedWidget, QMainWindow, QProgressBar
from PyQt5.QtCore import QTimer, Qt, QObject, pyqtSignal, QThread
from unyt import unyt_array, unyt_quantity
from scipy.spatial.transform import Rotation
from scipy.stats import iqr

try:
    import smplotlib
except:
    pass

#from lib.utils import create_sph_dataset, getBinNumber

Mbol_sun = 4.74 #4.83    

def create_sph_dataset(ad, pt, 
                       extra_fields = None, 
                       n_neighbours = 32, 
                       kernel = 'wendland6',
                       use_norm = True,
                       rotation_matrix = None
                      ):
    """Creates an sph field using the provided particle field, from the data contained
    in the ds dataset. The smoothing is performed by .add_sph_field. Both the number of
    neighbours and kernel can be selected. The returned field has name 'io'.

    Parameters
    ----------
    ad : yt.data_objects.selection_objects.region.YTRegion
        Original dataset.
    pt : str
        Particle type to be smoothened.
    extra_fields : list, optional
        Fields of the particles that will be added to the smoothened dataset.
        Default is None, which adds only mass, position and velocity of the particles.
        Any extra field passed, will be appended to this list.
    n_neighbours : int, optional
        Number of neighbour particles to include in the smoothing. Default is 32.
    kernel : str, optional
        Kernel used in the smoothing. Default is 'wendland6' (see yt Documentation).
    use_norm : bool, optional
        Value of use_sph_normalization. Default is True.

    Returns
    -------
    ds_sph : yt.dataset
        Dataset containing the sph field. The returned field has name 'io'.
    """
    if rotation_matrix is None:
        fields = ['particle_mass'] + [f'particle_position_{ax}' for ax in 'xyz'] + [f'particle_velocity_{ax}' for ax in 'xyz'] + ['particle_index']
        if extra_fields is not None:
            for e_field in extra_fields:
                field.append(e_field)
                
        data = {field: ad[pt, field] for field in fields}
    else:
        c = ad[pt,"Coordinates"] - unyt_array([512,512,512], 'kpc')
        v = ad[pt,"Velocities"]
        rotated_coords = unyt_array([np.dot(rotation_matrix, c[i,:]) for i in range(len(c))], ) + unyt_array([512,512,512], 'kpc')
        rotated_vels = unyt_array([np.dot(rotation_matrix, v[i,:]) for i in range(len(v))], )
        data = {
            'particle_mass' : ad[pt, 'Mass'],
            
            'particle_position_x' : rotated_coords[:,0],
            'particle_position_y' : rotated_coords[:,1],
            'particle_position_z' : rotated_coords[:,2],
            
            'particle_velocity_x' : rotated_vels[:,0],
            'particle_velocity_y' : rotated_vels[:,1],
            'particle_velocity_z' : rotated_vels[:,2],
            
            'particle_index' : ad[pt, 'particle_index']
        }

    ds_sph = yt.load_particles(data, data_source=ad)
    
    ds_sph.use_sph_normalization = use_norm
    
    ds_sph.add_sph_fields(n_neighbors=n_neighbours, kernel=kernel, sph_ptype='io')
    #ds_sph.add_deposited_particle_field(('io','particle_luminosity'), method="cic", kernel_name=kernel)
    return ds_sph
    
def getBinSize(data_on_axis):
    """Calculates optimum bin-length using the 
       Freedman-Diaconis rule (Freedman & Diaconis 1981).

       Parameters
       ----------
       data_on_axis : array
           Data that wants to be binned.
           Array of shape (n,)

       Returns
       -------
       h : float
           Optimum bin size to histogramize the data.
    """
    N = len(data_on_axis)
    h = 2*iqr(data_on_axis)/N**(1/3)
    return h

def getBinNumber(data_on_axis):
    """Calculates optimum bin-number based on the optimum bins-ize
       estimation by 'getBinSize'.

       Parameters
       ----------
       data_on_axis : array
           Data that wants to be binned.
           Array of shape (n,)

       Returns
       -------
       n : float
           Optimum bin number to histogramize the data.
    """
    n = int(np.rint((data_on_axis.max()-data_on_axis.min())/getBinSize(data_on_axis)))
    return n

def combinations_of_arrays(alpha, beta, gamma):
    # Generate all combinations of elements from alpha, beta, gamma
    import itertools
    combinations = list(itertools.product(alpha, beta, gamma))
    return combinations

def compute_surface_brightness(file_path, euler_angles, mtol, width, center):
    # Load data and create mock image
    ds = yt.load(file_path)
    ad = ds.all_data()
    t = float(ds.current_time.in_units('Gyr').value)
    
    rot_matrix = Rotation.from_euler(seq="xyz", angles=euler_angles, degrees=True).as_matrix()
    #self.rot_matrix = rot_matrix np.array([alpha, beta, gamma])
    
    
    ds_sph = create_sph_dataset(ad, "PartType4", n_neighbours=64, kernel="wendland2", rotation_matrix=rot_matrix) #64
    ad_sph = ds_sph.all_data()
    
    ML = unyt_quantity(mtol, 'Msun/Lsun')

    from yt import derived_field
    @derived_field(name="luminosity_density", sampling_type="particle", units="Lsun/pc**3")
    def _luminosity_density(field, data):
        try:
            return (
                data['io', 'density']/ML
            )
        except:
            return None
    
    ds_sph.add_field(("io", "luminosity"), function=_luminosity_density, sampling_type="particle", units="Lsun/pc**3")
    

    proj = yt.ProjectionPlot(ds_sph, [0,0,1], ('io','luminosity'), width=(width,'kpc'), center=unyt_array(center,'kpc'))
    proj.set_unit(('io','luminosity'), 'Lsun/pc**2')
    
    surflum = proj.to_fits_data().get_data('luminosity')
    

    surface_brightness = Mbol_sun + 21.572 - 2.5*np.log10((surflum).value)   #4.83
            
    return surflum, surface_brightness, rot_matrix, t













class AnalysisWorker(QThread):
    progress_updated = pyqtSignal(int)
    finished = pyqtSignal()

    def __init__(self, file_path, alpha, beta, gamma, mtol, width, center, distance, savefile = None):
        super().__init__()
        self.file_path = file_path
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.mtol = mtol
        self.width = width
        self.center = center
        self.distance = distance
        self.savefile = savefile

        print(self.file_path,
        self.alpha,
        self.beta,
        self.gamma,
        self.mtol,
        self.width,
        self.center,
        self.savefile)
    
    def run(self):
        if self.file_path.endswith(".hdf5"):
            self.run_on_hdf5()
        else:
            self.run_on_folder()
    
    def run_on_folder(self):
        euler_combinations = combinations_of_arrays(self.alpha, self.beta, self.gamma)
        files = np.array([file for file in os.listdir(self.file_path) if file.endswith('.hdf5')])
        files.sort()
        total_files = len(files)

        self.sb_max = {str(key): np.zeros_like(files, dtype="float") for key in euler_combinations}
        self.sb_max['snapshot'] = files
        self.sb_max['time'] = np.zeros_like(files, dtype="float")
        
        self.observed_area = {str(key): np.zeros_like(files, dtype="float") for key in euler_combinations}
        self.observed_area['snapshot'] = files
        self.observed_area['time'] = np.zeros_like(files, dtype="float")

     
        
        for i, file in enumerate(files):
            for j, angles in enumerate(euler_combinations):
                sf, sb, _, t = compute_surface_brightness(self.file_path+"/"+file, angles, self.mtol, self.width, self.center)

                s = (unyt_quantity(self.width, 'rad*kpc')/unyt_quantity(self.distance,'Mpc')).in_units('arcsec')
                l = unyt_quantity(10,'arcsec')
                tot_pix, _ = sb.shape
                npix = int(np.ceil(tot_pix * l/s))
                nside, res = npix//2, npix%2

                indexes = (ii, jj) = np.unravel_index(sb.argmin(), sb.shape)

                sf_mean = sf[ii-nside:ii+nside+res,jj-nside:jj+nside+res].mean()
                sb_mean =  Mbol_sun + 21.572 - 2.5*np.log10(sf_mean.value)  #4.83
                
                self.sb_max['time'][i] = t
                self.observed_area['time'][i] = t
                
                self.sb_max[str(angles)][i] = sb_mean
                self.observed_area[str(angles)][i] = len(sb[sb <= 30])*(s/tot_pix).value**2
               

            
            progress_percentage = int((i + 1) * 100 / len(files))
            self.progress_updated.emit(progress_percentage)

                
        pd_sb = pd.DataFrame(self.sb_max)
        pd_sb.to_csv(f'{self.savefile}_sbmax.csv', index=False)
        
        pd_oa = pd.DataFrame(self.observed_area)
        pd_oa.to_csv(f'{self.savefile}_oa.csv', index=False)
        
        self.finished.emit()

    
    def run_on_hdf5(self):
        euler_combinations = combinations_of_arrays(self.alpha, self.beta, self.gamma)
        total_combinations = len(euler_combinations)
        self.derived_quantities = {'euler_angles' : np.empty((total_combinations,), dtype="<U100"),
                                   'sb_max' : np.empty((total_combinations,), dtype="float"),
                                   'obs_area' : np.empty((total_combinations,), dtype="float")
                                  }
        
        for i, angles in enumerate(euler_combinations):
            sf, sb, _, _ = compute_surface_brightness(self.file_path, angles, self.mtol, self.width, self.center)

            s = (unyt_quantity(self.width, 'rad*kpc')/unyt_quantity(self.distance,'Mpc')).in_units('arcsec')
            l = unyt_quantity(10,'arcsec')
            tot_pix, _ = sb.shape
            npix = int(np.ceil(tot_pix * l/s))
            nside, res = npix//2, npix%2

            indexes = (ii, jj) = np.unravel_index(sb.argmin(), sb.shape)

            sf_mean = sf[ii-nside:ii+nside+res,jj-nside:jj+nside+res].mean()
            sb_mean =  Mbol_sun + 21.572 - 2.5*np.log10(sf_mean.value)  #4.83
            
            observable_pixels = len(sb[sb <= 30])*(s/tot_pix).value**2
            
            self.derived_quantities['euler_angles'][i] = str(angles)
            self.derived_quantities['sb_max'][i] = sb_mean
            self.derived_quantities['obs_area'][i] = observable_pixels
          
            
            progress_percentage = int((i + 1) * 100 / total_combinations)
            self.progress_updated.emit(progress_percentage)

        pd_dq = pd.DataFrame(self.derived_quantities)
        pd_dq.to_csv('data.csv', index=False)
        self.finished.emit()



class AnalysisWindow(QWidget):
    def __init__(self, file_path, center, width, distance, mtol):
        super().__init__()
        self.file_path = file_path
        self.center = center
        self.width = width
        self.distance = distance
        self.mtol = mtol

        self.setWindowTitle('Analysis')
        self.setGeometry(100, 100, 600, 300)
        
        self.alpha_label = QLabel('Alpha (degrees):')
        self.alpha_edit = QLineEdit()

        self.beta_label = QLabel('Beta (degrees):')
        self.beta_edit = QLineEdit()

        self.gamma_label = QLabel('Gamma (degrees):')
        self.gamma_edit = QLineEdit()

        self.mtol_label = QLabel('ML-Ratio:')
        self.mtol_edit = QLineEdit()

    
        self.folder_label = QLabel('Folder/File Path:')
        self.folder_edit = QLineEdit()
        self.browse_button = QPushButton('Browse')
        self.browse_button.clicked.connect(self.browseFolder)
        
        self.analysis_button = QPushButton('Perform Analysis')
        self.analysis_button.clicked.connect(self.performAnalysis)

        self.progress_bar = QProgressBar()
        self.progress_bar.setTextVisible(True)
        self.progress_bar.setStyleSheet("""
            QProgressBar {
                border: 1px solid black;
                border-radius: 8px;
                text-align: center;
                height: 10px;  /* Adjust the height of the progress bar */
            }
            QProgressBar::chunk {
                background-color: #37c9e1;
                width: 10px;  /* Adjust the width here to make the progress bar thinner */
            }
        """)

        main_layout = QVBoxLayout()
        
        file_layout = QHBoxLayout()
        file_layout.addWidget(self.folder_label)
        file_layout.addWidget(self.folder_edit)
        file_layout.addWidget(self.browse_button)
        
        others_layout = QVBoxLayout()
        others_layout.addWidget(self.alpha_label)
        others_layout.addWidget(self.alpha_edit)
        others_layout.addWidget(self.beta_label)
        others_layout.addWidget(self.beta_edit)
        others_layout.addWidget(self.gamma_label)
        others_layout.addWidget(self.gamma_edit)
        others_layout.addWidget(self.mtol_label)
        others_layout.addWidget(self.mtol_edit)
        others_layout.addWidget(self.progress_bar)
        others_layout.addWidget(self.analysis_button)

        main_layout.addLayout(file_layout)
        main_layout.addLayout(others_layout)
        
        self.setLayout(main_layout)

    def performAnalysis(self):
        alpha = [float(x) for x in self.alpha_edit.text().split(',')]
        beta = [float(x) for x in self.beta_edit.text().split(',')]
        gamma = [float(x) for x in self.gamma_edit.text().split(',')]

        file_path = self.folder_edit.text() if self.folder_edit.text() else self.file_path
        mtol = float(self.mtol_edit.text()) if self.mtol_edit.text() else self.mtol

        savefile = None if file_path.endswith('hdf5') else file_path.split("/")[-2]
        
        self.progress_bar.setValue(0)
        self.analysis_worker = AnalysisWorker(file_path, alpha, beta, gamma, mtol, self.width, self.center, self.distance, savefile)
        self.analysis_worker.progress_updated.connect(self.updateProgress)
        self.analysis_worker.finished.connect(self.onAnalysisFinished)
        self.analysis_worker.start()

    def updateProgress(self, percentage):
        self.progress_bar.setValue(percentage)
    
    def onAnalysisFinished(self):
        print("Analysis completed.")
        self.progress_bar.setValue(100)
        # Perform any post-analysis actions here

    def browseFolder(self):
        folder_path = QFileDialog.getExistingDirectory(self, 'Select Folder')
        if folder_path:
            self.folder_edit.setText(folder_path)












class MockImageApp(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()
        
    def initUI(self):
        self.surface_brightness = None
        self.file_path = None
        
        self.setWindowTitle('Mock Image Options')
        self.setGeometry(100, 100, 600, 300)


        
        # Create input fields for each argument
        self.file_label = QLabel('File Path:')
        self.file_edit = QLineEdit()
        self.file_button = QPushButton('Browse')
        self.file_button.clicked.connect(self.browseFile)
        
        self.alpha_label = QLabel('Alpha (degrees):')
        self.alpha_edit = QLineEdit()
        
        self.beta_label = QLabel('Beta (degrees):')
        self.beta_edit = QLineEdit()
        
        self.gamma_label = QLabel('Gamma (degrees):')
        self.gamma_edit = QLineEdit()
        
        self.distance_label = QLabel('Distance (Mpc):')
        self.distance_edit = QLineEdit()
        
        self.center_label = QLabel('Center (x, y, z) in kpc:')
        self.center_edit = QLineEdit()
        
        self.width_label = QLabel('Width (kpc):')
        self.width_edit = QLineEdit()
        
        self.contours_label = QLabel(r'Contours (mag/arcsec^2):')
        self.contours_edit = QLineEdit()
        
        self.clims_label = QLabel('Color Limits (clow, chigh):')
        self.clims_edit = QLineEdit()
        
        self.mtol_label = QLabel('ML-Ratio:')
        self.mtol_edit = QLineEdit()
        
        self.create_button = QPushButton('Create Mock Image')
        self.create_button.clicked.connect(self.createMockImage)
        self.create_button.setFixedSize(280, 30)  # Set fixed size for button
        
        self.to_fits_button = QPushButton('To FITS')
        self.to_fits_button.clicked.connect(self.saveToFITS)
        self.to_fits_button.setFixedSize(280, 30)  # Set fixed size for button
        
        self.fits_name_label = QLabel('FITS File:')
        self.fits_name_edit = QLineEdit()

        # Checkbox for histogram option
        self.histogram_checkbox = QCheckBox('Include Histogram')
        self.histogram_checkbox.stateChanged.connect(self.toggleHistogram)
        self.show_histogram = False
        
        # Analysis tab
        self.analysis_button = QCheckBox('Open Analysis')
        self.analysis_button.clicked.connect(self.openAnalysisWindow)
        
        # Arrange widgets in vertical layout
        main_layout = QVBoxLayout()
        
        # File Path line
        file_layout = QHBoxLayout()
        file_layout.addWidget(self.file_label)
        file_layout.addWidget(self.file_edit)
        file_layout.addWidget(self.file_button)
        
        # Left column (Alpha, Beta, Gamma, Distance)
        left_layout = QVBoxLayout()
        left_layout.addWidget(self.alpha_label)
        left_layout.addWidget(self.alpha_edit)
        left_layout.addWidget(self.beta_label)
        left_layout.addWidget(self.beta_edit)
        left_layout.addWidget(self.gamma_label)
        left_layout.addWidget(self.gamma_edit)
        left_layout.addWidget(self.distance_label)
        left_layout.addWidget(self.distance_edit)
        
        # Right column (Center, Width, Contours, Color Limits)
        right_layout = QVBoxLayout()
        right_layout.addWidget(self.center_label)
        right_layout.addWidget(self.center_edit)
        right_layout.addWidget(self.width_label)
        right_layout.addWidget(self.width_edit)
        right_layout.addWidget(self.contours_label)
        right_layout.addWidget(self.contours_edit)
        right_layout.addWidget(self.clims_label)
        right_layout.addWidget(self.clims_edit)

        # Mass-to-Light Ratio field
        mtol_layout = QHBoxLayout()
        mtol_layout.addWidget(self.mtol_label)
        mtol_layout.addWidget(self.mtol_edit)
        
        # FITS file name layout
        fits_name_layout = QHBoxLayout()
        fits_name_layout.addWidget(self.fits_name_label)
        fits_name_layout.addWidget(self.fits_name_edit)
        
        # Create Mock Image and To FITS buttons layout
        buttons_layout = QHBoxLayout()
        buttons_layout.addWidget(self.create_button)
        buttons_layout.addWidget(self.to_fits_button)

        # Add histogram checkbox to the layout
        checkbox_layout = QHBoxLayout()
        checkbox_layout.addWidget(self.histogram_checkbox)
        checkbox_layout.addWidget(self.analysis_button)

        
        # Columns layout
        columns_layout = QHBoxLayout()
        columns_layout.addLayout(left_layout)
        columns_layout.addLayout(right_layout)


        # Add widgets to main layout
        main_layout.addLayout(file_layout)
        main_layout.addLayout(columns_layout)
        main_layout.addLayout(mtol_layout)
        main_layout.addLayout(fits_name_layout)
        main_layout.addLayout(checkbox_layout)
        main_layout.addLayout(buttons_layout)
        
        self.setLayout(main_layout)
        self.show()

        self.fade_timer = QTimer()
        self.fade_timer.setInterval(50)  # Timer interval (50 ms for smooth transition)
        self.fade_timer.timeout.connect(self.fadeOutButtons)

    def toggleHistogram(self, state):
        self.show_histogram = state == Qt.Checked
        pass

    def openAnalysisWindow(self):
        if self.analysis_button.isChecked():  # Check if the analysis checkbox is checked
            self.file_path = self.file_edit.text() if self.file_edit.text() else None
            self.center = [float(x) for x in self.center_edit.text().split(',')] if self.center_edit.text() else [512,512,512]
            self.width = float(self.width_edit.text()) if self.width_edit.text() else 200
            self.distance = float(self.distance_edit.text()) if self.distance_edit.text() else 45
            self.mtol = float(self.mtol_edit.text()) if self.mtol_edit.text() else None
            
            self.analysis_window = AnalysisWindow(self.file_path, self.center, self.width, self.distance, self.mtol)  # Create an instance of AnalysisWindow
            self.analysis_window.show()  # Show the AnalysisWindow
        else:
            if hasattr(self, 'analysis_window'):  # Check if the analysis window instance exists
                self.analysis_window.close() 


    def fadeOutButtons(self):
        # Restore the original button styles after the fade-out effect
        self.create_button.setStyleSheet("")
        self.to_fits_button.setStyleSheet("")
        self.fade_timer.stop()
        pass
        
    def browseFile(self):
        options = QFileDialog.Options()
        file_path, _ = QFileDialog.getOpenFileName(self, "Select File", "", "All Files (*);;Python Files (*.py)", options=options)
        if file_path:
            self.file_edit.setText(file_path)
        pass

    def plotImage(self):
        from mpl_toolkits.axes_grid1 import AxesGrid
        fig = plt.figure()
        grid = AxesGrid(
            fig,
            (0.075, 0.075, 0.85, 0.85),
            nrows_ncols=(1, 1),
            axes_pad=0.05,
            label_mode="L",
            share_all=True,
            cbar_location="right",
            cbar_mode="single",
            cbar_size="5%",
            cbar_pad="0%",
        )
    
        
        s = (unyt_quantity(self.width/2, 'rad*kpc')/unyt_quantity(self.distance,'Mpc')).in_units('arcmin')
        n,m = self.surface_brightness.shape
        x = np.linspace(-s+s/n, s-s/n,n)
        X, Y = np.meshgrid(x,x)
    
        
        ax = grid[0]
        
        ax.set_title(f"t={self.t:.1f} Gyr Mock image with M/L = {self.mtol}\nsb_max={self.surface_brightness.min():.2f} mag·arcsec^-2\n")
        
        ax.set_xlabel(r"$\alpha$ [arcmin]")
        ax.set_ylabel(r"$\delta$ [arcmin]")
        
        cmyt.arbre_r.set_bad(cmyt.arbre_r.get_over())
        
        im = ax.pcolormesh(X,Y,self.surface_brightness, cmap=cmyt.arbre_r)
        im.set_clim(self.chigh, self.clow)
        
        cbar = grid.cbar_axes[0].colorbar(im)
        
        ticks = np.arange(self.chigh, self.clow, 2, dtype="int")
        tick_labels = [f'{tick:.0f}' for tick in ticks] 
        cbar.set_ticks(ticks)
        cbar.set_ticklabels(tick_labels)
        cbar.set_label('Surface Brightness [mag/arcsec^2]')
        
        if self.contour_levels is not None:    
            contour = ax.contour(X, Y, self.surface_brightness, levels=self.contour_levels, colors='red',linewidths=0.7)
            for level in self.contour_levels:
                cbar.ax.axhline(level, color='red', linewidth=1.5, linestyle='-')
        
        plt.show()
        #return grid

    def plotHistogram(self):
        fig, ax = plt.subplots(figsize=(5.1,5.1))

        plt.title(f"Pixel Histogram")
        
        ax.set_xlabel(f"Surface Brightness [mag/arcsec^2]")
        ax.set_ylabel(f"Counts")
        
        mask = (~np.isinf(self.surface_brightness))
        ax.hist(self.surface_brightness[mask], log=True, bins=getBinNumber(self.surface_brightness[mask]))

        quantiles = np.quantile(self.surface_brightness[mask], [0.05, 0.5, 0.95])
        ax.axvline(x = quantiles[0], ymin=0, ymax=1E5, color="blue", linewidth=0.5)
        ax.axvline(x = quantiles[1], ymin=0, ymax=1E5, color="blue", linewidth=0.5)
        ax.axvline(x = quantiles[2], ymin=0, ymax=1E5, color="blue", linewidth=0.5)

        if self.contour_levels is not None:    
            for level in self.contour_levels:
                ax.axvline(x = level, ymin=0, ymax=1E5, color="red", linewidth=0.5)

        plt.show()
        #return grid

    def createMockImage(self):
        success = True
        try:
            # Get input values from QLineEdit fields
            self.file_path = self.file_edit.text()
            alpha = float(self.alpha_edit.text())
            beta = float(self.beta_edit.text())
            gamma = float(self.gamma_edit.text())
            self.center = [float(x) for x in self.center_edit.text().split(',')] if self.center_edit.text() else [512,512,512]
            self.width = float(self.width_edit.text()) if self.width_edit.text() else 200
            self.distance = float(self.distance_edit.text()) 
            self.mtol = float(self.mtol_edit.text())
            self.contour_levels = [float(x) for x in self.contours_edit.text().split(',')] if self.contours_edit.text() else None
            clims = [float(x) for x in self.clims_edit.text().split(',')] if self.clims_edit.text() else [29, 40]
            self.chigh, self.clow = clims[0], clims[1]
    
            _, surface_brightness, rot_matrix, t = compute_surface_brightness(self.file_path, (alpha,beta,gamma), self.mtol, self.width, self.center)
            self.surface_brightness = surface_brightness
            self.rot_matrix = rot_matrix
            self.t = t

            if self.show_histogram is True:
                self.plotImage()
                self.plotHistogram()
            else:
                self.plotImage()
            
            
        except:
            success = False

        if success:
            self.create_button.setStyleSheet("background-color: green; color: white;")
        else:
            self.create_button.setStyleSheet("background-color: red; color: white;")

        self.fade_timer.start(5000)
        pass

    def saveToFITS(self):
        success = True
        try:
            fits_name = self.fits_name_edit.text().strip()
            if not fits_name:
                fits_name = "mock_image.default.fits"  # Default FITS file name if not specified
            
            if self.surface_brightness is not None:
                import datetime
                from astropy.io import fits
                hdr = fits.Header()
                hdr['SIMPLE'] = 'T'
                hdr['BITPIX'] = -64  # 64-bit floating point
                hdr['NAXIS'] = 2
                hdr['NAXIS1'] = self.surface_brightness.shape[1]
                hdr['NAXIS2'] = self.surface_brightness.shape[0]
                hdr['OBJECT'] = 'N-BODY SIMULATION MOCK IMAGE'   # Name of the observed object
                hdr['CODE'] = 'SWIFT'
                current_date = datetime.date.today()
                hdr['DATE-OBS'] = current_date.strftime('%Y-%m-%d')  # Date of observation (YYYY-MM-DD)
                hdr['AUTHOR'] = 'Asier Lambarri Martinez'   # Name of the observer
                hdr['BUNIT'] = 'mag/arcsec^2'            # Unit of the pixel values (Analog-to-Digital Units)
                hdr['FILTER'] = 'BOLOMETRIC'
                hdr['PHYSICAL_IMAGE_SIZE'] = f"{float(self.width_edit.text())} kpc" if self.width_edit.text() else f"512 kpc"             
                hdr['DISTANCE_TO_OBJECT'] = f"{float(self.distance_edit.text())} Mpc"             
                p = np.dot(self.rot_matrix, [0,0,1])
                hdr['PROJECTION'] = ', '.join(map(str, p))
                hdr['SNAPSHOT_TIME'] = f"{self.t} Gyr"
                hdr['COMMENT'] = ''
                hdu = fits.PrimaryHDU(self.surface_brightness, header=hdr)
                
                hdu.writeto(fits_name, overwrite=True)
    
            else:
                print("No surface brightness data available to save.")
                success = False
        except:
            success = False
        if success:
            self.to_fits_button.setStyleSheet("background-color: green; color: white;")
        else:
            self.to_fits_button.setStyleSheet("background-color: red; color: white;")
            
        self.fade_timer.start(5000)
        pass













if __name__ == '__main__':
    app = QApplication(sys.argv)
    mock_image_app = MockImageApp()
    sys.exit(app.exec_())

