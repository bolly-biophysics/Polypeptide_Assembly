#!/usr/bin/env python3
"""
WT-MTD Free Energy Surface Reconstruction Tool (Unit conversion: kJ/mol → kcal/mol)
Modification: Added contour line and EPS output support
New features:
1. Coordinate axis range restriction
2. Color bar upper limit setting
3. Mark initial coordinate points
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import os

class WT_FES_Reconstructor:
    def __init__(self, gamma, input_unit='kJ/mol', output_unit='kcal/mol', 
                 kT=None, temp=300.0):
        """
        Initialize reconstructor
        
        Parameters:
        gamma: bias factor
        input_unit: Unit of input negbias file ('kJ/mol' or 'kcal/mol')
        output_unit: Unit of output free energy ('kJ/mol' or 'kcal/mol')
        kT: k_B*T value (energy units), overrides temp if provided
        temp: Temperature (K)
        """
        self.gamma = float(gamma)
        self.input_unit = input_unit
        self.output_unit = output_unit
        
        # Calculate conversion factor
        self.conversion_factor = self.gamma / (self.gamma - 1.0)
        print(f"WT-MTD conversion factor γ/(γ-1) = {self.conversion_factor:.6f}")
        
        # Unit conversion factor
        if input_unit == output_unit:
            self.unit_conversion_factor = 1.0
        elif input_unit == 'kJ/mol' and output_unit == 'kcal/mol':
            self.unit_conversion_factor = 0.239006  # 1 kJ/mol = 0.239006 kcal/mol
        elif input_unit == 'kcal/mol' and output_unit == 'kJ/mol':
            self.unit_conversion_factor = 4.184  # 1 kcal/mol = 4.184 kJ/mol
        else:
            raise ValueError(f"Unsupported unit conversion: {input_unit} → {output_unit}")
        
        print(f"Unit conversion factor: 1 {input_unit} = {self.unit_conversion_factor} {output_unit}")
        
        # Calculate kT value (for information display)
        if kT is not None:
            self.kT_input = float(kT)  # User-provided kT value (input units)
            # Convert to output units
            self.kT_output = self.kT_input * self.unit_conversion_factor
            print(f"Using provided kT value: {self.kT_input:.4f} {input_unit} = {self.kT_output:.4f} {output_unit}")
        else:
            # Calculate kT based on input units
            if input_unit == 'kcal/mol':
                kB = 0.0019872041  # kcal/(mol·K)
            else:  # kJ/mol
                kB = 0.008314462618  # kJ/(mol·K)
            self.kT_input = kB * temp
            self.kT_output = self.kT_input * self.unit_conversion_factor
            print(f"kT calculated from temperature {temp}K: {self.kT_input:.4f} {input_unit} = {self.kT_output:.4f} {output_unit}")
        
    def reconstruct_from_negbias(self, negbias_file, output_file=None, 
                                 shift_min_to_zero=True):
        """
        Reconstruct free energy surface from negbias file and perform unit conversion
        """
        print(f"\n{'='*60}")
        print(f"Step 1: Read negbias file: {negbias_file}")
        print(f"Input unit: {self.input_unit}, Output unit: {self.output_unit}")
        
        # Read data
        data = np.loadtxt(negbias_file)
        
        # Check data dimensions
        if data.shape[1] < 3:
            raise ValueError(f"Data requires at least 3 columns, but has {data.shape[1]} columns")
        
        # Extract data (values in negbias file are negative bias potential, unit is input_unit)
        cv1 = data[:, 0]
        cv2 = data[:, 1]
        negbias_input = data[:, 2]  # Unit: input_unit
        
        print(f"Number of data points: {len(negbias_input)}")
        print(f"negbias range (input units): [{np.min(negbias_input):.4f}, {np.max(negbias_input):.4f}] {self.input_unit}")
        
        # Step 1: Apply WT-MTD conversion factor to obtain correct free energy (still in input units)
        print(f"\nStep 2: Apply WT-MTD conversion factor {self.conversion_factor:.6f}")
        fes_input_units = self.conversion_factor * negbias_input
        
        print(f"Free energy range after conversion (input units): [{np.min(fes_input_units):.4f}, {np.max(fes_input_units):.4f}] {self.input_unit}")
        
        # Step 2: Unit conversion
        print(f"\nStep 3: Unit conversion: {self.input_unit} → {self.output_unit}")
        print(f"Unit conversion factor: {self.unit_conversion_factor}")
        fes_output_units = fes_input_units * self.unit_conversion_factor
        
        print(f"Free energy range after conversion (output units): [{np.min(fes_output_units):.4f}, {np.max(fes_output_units):.4f}] {self.output_unit}")
        
        # Step 3: Optional: shift minimum to zero
        if shift_min_to_zero:
            fes_min = np.min(fes_output_units)
            fes_final = fes_output_units - fes_min
            print(f"Shift minimum to zero (original minimum: {fes_min:.4f} {self.output_unit})")
        else:
            fes_final = fes_output_units
        
        # Prepare output filename
        if output_file is None:
            base_name = os.path.splitext(negbias_file)[0]
            unit_suffix = self.output_unit.replace('/', '_')
            output_file = f"{base_name}_fes_{unit_suffix}.dat"
        
        # Create output data
        output_data = np.column_stack([cv1, cv2, fes_final])
        
        # Generate detailed file header
        header = f"""# Free Energy Surface reconstructed from WT-MTD negbias
# Source file: {os.path.abspath(negbias_file)}
# 
# RECONSTRUCTION PARAMETERS:
#   Bias factor (gamma) = {self.gamma}
#   WT-MTD conversion factor = γ/(γ-1) = {self.conversion_factor:.8f}
#   
# UNIT CONVERSION:
#   Input negbias unit = {self.input_unit}
#   Output FES unit = {self.output_unit}
#   Unit conversion factor = {self.unit_conversion_factor:.8f}
#   kT (input units) = {self.kT_input:.6f} {self.input_unit}
#   kT (output units) = {self.kT_output:.6f} {self.output_unit}
#   
# POST-PROCESSING:
#   Minimum shifted to zero = {'Yes' if shift_min_to_zero else 'No'}
#   
# COLUMNS:
#   1. CV1
#   2. CV2  
#   3. Free_Energy({self.output_unit})
#
# CONVERSION SUMMARY:
#   FES_output = [negbias_input × γ/(γ-1)] × (unit_conversion_factor)
#               = [negbias_input × {self.conversion_factor:.6f}] × {self.unit_conversion_factor:.6f}
"""
        
        # Save file
        np.savetxt(output_file, output_data, header=header, fmt='%14.8f')
        
        print(f"\nStep 4: Save result to {output_file}")
        print(f"Final free energy range: [{np.min(fes_final):.4f}, {np.max(fes_final):.4f}] {self.output_unit}")
        
        # Analyze minimum
        min_idx = np.argmin(fes_final)
        print(f"\nGlobal minimum:")
        print(f"  CV1 = {cv1[min_idx]:.6f}")
        print(f"  CV2 = {cv2[min_idx]:.6f}")
        print(f"  FES = {fes_final[min_idx]:.6f} {self.output_unit}")
        
        # If needed, output values in original units as well
        if shift_min_to_zero:
            fes_input_min = np.min(fes_input_units)
            fes_input_zeroed = fes_input_units - fes_input_min
            print(f"  FES (original units, zero-shifted) = {fes_input_zeroed[min_idx]:.6f} {self.input_unit}")
        
        return cv1, cv2, fes_final, output_file
    
    def plot_fes_heatmap(self, cv1, cv2, fes, plot_file=None, 
                         cmap='RdYlBu_r', dpi=300, 
                         add_min_marker=True, add_grid=False,
                         min_label_format='cv_only',
                         add_contour=True, contour_levels=10,
                         contour_colors='black', contour_alpha=0.7,
                         contour_linewidths=0.8, output_format='png',
                         # New parameters
                         cv1_range=None, cv2_range=None,
                         vmax=None, initial_coords=None):
        """
        Plot free energy surface heatmap, optionally with contour lines
        
        Parameters:
        add_min_marker: Whether to mark the minimum point
        add_grid: Whether to add grid
        min_label_format: Minimum label format
            - 'cv_only': Show only CV coordinates
            - 'fes_only': Show only free energy value
            - 'both': Show CV coordinates and free energy value
            - 'cv_fes': Show "CV=(x,y), F=value"
        add_contour: Whether to add contour lines
        contour_levels: Number of contour levels or list of specific values
        contour_colors: Contour line color
        contour_alpha: Contour line transparency
        contour_linewidths: Contour line width
        output_format: Output image format ('png', 'eps', 'pdf', 'svg')
        
        New parameters:
        cv1_range: CV1 display range [min, max]
        cv2_range: CV2 display range [min, max]
        vmax: Color bar maximum value (free energy upper limit)
        initial_coords: List of initial coordinate points [(cv1, cv2), ...] for marking starting points
        """
        print(f"\n{'='*60}")
        print("Step 5: Plot free energy surface heatmap")
        print(f"Minimum label format: {min_label_format}")
        print(f"Add contour lines: {'Yes' if add_contour else 'No'}")
        print(f"Contour levels: {contour_levels}")
        print(f"Output format: {output_format}")
        
        # Print new parameter information
        if cv1_range is not None:
            print(f"CV1 display range: [{cv1_range[0]}, {cv1_range[1]}]")
        if cv2_range is not None:
            print(f"CV2 display range: [{cv2_range[0]}, {cv2_range[1]}]")
        if vmax is not None:
            print(f"Color bar maximum: {vmax} {self.output_unit}")
        if initial_coords is not None:
            print(f"Number of initial coordinate points to mark: {len(initial_coords)}")
            for i, coord in enumerate(initial_coords):
                print(f"  Point {i+1}: ({coord[0]}, {coord[1]})")
        
        # Reshape to grid (assuming data is regular)
        cv1_vals = np.unique(cv1)
        cv2_vals = np.unique(cv2)
        
        if len(cv1_vals) * len(cv2_vals) != len(fes):
            print("⚠️ Data point count mismatch, may not be a regular grid")
            print(f"Unique CV1 values: {len(cv1_vals)}, Unique CV2 values: {len(cv2_vals)}")
            print(f"Total data points: {len(fes)}")
            print("Attempting to reshape grid...")
        
        try:
            # Attempt to reshape
            cv1_grid = cv1.reshape(len(cv2_vals), len(cv1_vals))
            cv2_grid = cv2.reshape(len(cv2_vals), len(cv1_vals))
            fes_grid = fes.reshape(len(cv2_vals), len(cv1_vals))
            
            # Find minimum position
            min_idx_2d = np.unravel_index(np.argmin(fes_grid), fes_grid.shape)
            min_cv1 = cv1_grid[min_idx_2d]
            min_cv2 = cv2_grid[min_idx_2d]
            min_fes = fes_grid[min_idx_2d]
            
            # Generate label based on format
            if min_label_format == 'cv_only':
                min_label = f'Min: ({min_cv1:.2f}, {min_cv2:.2f})'
            elif min_label_format == 'fes_only':
                min_label = f'Min F: {min_fes:.2f} {self.output_unit}'
            elif min_label_format == 'both':
                min_label = f'Min: ({min_cv1:.2f}, {min_cv2:.2f}), F: {min_fes:.2f} {self.output_unit}'
            elif min_label_format == 'cv_fes':
                min_label = f'CV=({min_cv1:.2f}, {min_cv2:.2f}), F={min_fes:.2f} {self.output_unit}'
            else:
                min_label = f'Min: ({min_cv1:.2f}, {min_cv2:.2f})'
            
            # Create figure
            fig, ax = plt.subplots(figsize=(10, 8))
            
            # Determine display range
            if cv1_range is not None:
                xmin, xmax = cv1_range
            else:
                xmin, xmax = np.min(cv1_grid), np.max(cv1_grid)
                
            if cv2_range is not None:
                ymin, ymax = cv2_range
            else:
                ymin, ymax = np.min(cv2_grid), np.max(cv2_grid)
            
            extent = [xmin, xmax, ymin, ymax]
            
            # Determine color bar range
            if vmax is not None:
                vmin = np.min(fes_grid)
                levels = np.linspace(vmin, vmax, 100)  # Use specified vmax
            else:
                vmin = np.min(fes_grid)
                vmax = np.max(fes_grid)
                levels = 100  # Use many levels for smooth gradient
            
            # Use contourf to draw heatmap
            contourf = ax.contourf(cv1_grid, cv2_grid, fes_grid, 
                                  levels=levels, cmap=cmap, 
                                  extent=extent, origin='lower',
                                  vmin=vmin, vmax=vmax)
            
            # Add contour lines
            if add_contour:
                if isinstance(contour_levels, int):
                    # Automatically generate contour levels
                    contour_plot = ax.contour(cv1_grid, cv2_grid, fes_grid, 
                                             levels=contour_levels, 
                                             colors=contour_colors,
                                             alpha=contour_alpha,
                                             linewidths=contour_linewidths)
                    # Add contour labels
                    ax.clabel(contour_plot, inline=True, fontsize=12, fmt='%.1f')
                else:
                    # Use specified contour levels
                    contour_plot = ax.contour(cv1_grid, cv2_grid, fes_grid, 
                                             levels=contour_levels, 
                                             colors=contour_colors,
                                             alpha=contour_alpha,
                                             linewidths=contour_linewidths)
                    ax.clabel(contour_plot, inline=True, fontsize=11, fmt='%.1f')
            
            # Set axis ranges
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(ymin, ymax)
            
            # Set axis labels
            ax.set_xlabel('CV1', fontsize=14)
            ax.set_ylabel('CV2', fontsize=14)
            ax.set_title(f'Free Energy Surface (γ={self.gamma}, unit: {self.output_unit})', 
                        fontsize=16, pad=20)
            
            # Add color bar
            cbar = plt.colorbar(contourf, ax=ax, pad=0.03)
            cbar.set_label(f'Free Energy ({self.output_unit})', fontsize=14)
            
            # Mark minimum
            if add_min_marker:
                ax.plot(min_cv1, min_cv2, 'r*', markersize=15, 
                       markeredgecolor='black', markeredgewidth=1,
                       label=min_label)
                ax.legend(loc='upper right', fontsize=11, framealpha=0.9)
            
            # Mark initial coordinate points
            if initial_coords is not None:
                for i, (init_cv1, init_cv2) in enumerate(initial_coords):
                    # Check if coordinate point is within display range
                    if xmin <= init_cv1 <= xmax and ymin <= init_cv2 <= ymax:
                        # Use different marker styles
                        ax.plot(init_cv1, init_cv2, 'go', markersize=10,
                               markeredgecolor='black', markeredgewidth=1,
                               alpha=0.8, label=f'Initial point {i+1}' if i == 0 else "")
                    else:
                        print(f"Warning: Initial coordinate point ({init_cv1}, {init_cv2}) is outside display range")
                
                # If points were marked, add legend
                if len(initial_coords) > 0:
                    # Get existing legend handles and labels
                    handles, labels = ax.get_legend_handles_labels()
                    # Add new legend
                    ax.legend(handles, labels, loc='upper right', fontsize=11, framealpha=0.9)
            
            # Add grid (optional)
            if add_grid:
                ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
            
            plt.tight_layout()
            
            # Generate output filename
            if plot_file is None:
                unit_suffix = self.output_unit.replace('/', '_')
                plot_file = f'fes_heatmap_gamma{self.gamma}_{unit_suffix}'
            else:
                # Remove existing extension, add new extension
                plot_file = os.path.splitext(plot_file)[0]
            
            # Save image according to format
            output_formats = output_format.lower().split(',')
            saved_files = []
            
            for fmt in output_formats:
                fmt = fmt.strip()
                if fmt not in ['png', 'eps', 'pdf', 'svg', 'jpg', 'tiff']:
                    print(f"Warning: Unsupported format '{fmt}', skipping")
                    continue
                
                output_filename = f"{plot_file}.{fmt}"
                if fmt == 'eps':
                    plt.savefig(output_filename, dpi=dpi, bbox_inches='tight', 
                               facecolor='white', format='eps')
                elif fmt == 'pdf':
                    plt.savefig(output_filename, dpi=dpi, bbox_inches='tight', 
                               facecolor='white', format='pdf')
                else:
                    plt.savefig(output_filename, dpi=dpi, bbox_inches='tight', 
                               facecolor='white', format=fmt)
                
                saved_files.append(output_filename)
                print(f"Image saved as: {output_filename} ({fmt.upper()})")
            
            plt.show()
            
            return fig, ax, saved_files
            
        except Exception as e:
            print(f"Heatmap plotting failed: {e}")
            print("Attempting scatter heatmap...")
            
            # Fallback: scatter heatmap
            fig, ax = plt.subplots(figsize=(13, 3))
            
            # Filter data points (based on display range)
            mask = np.ones_like(cv1, dtype=bool)
            if cv1_range is not None:
                mask &= (cv1 >= cv1_range[0]) & (cv1 <= cv1_range[1])
            if cv2_range is not None:
                mask &= (cv2 >= cv2_range[0]) & (cv2 <= cv2_range[1])
            
            cv1_filtered = cv1[mask]
            cv2_filtered = cv2[mask]
            fes_filtered = fes[mask]
            
            # Determine color bar range
            if vmax is not None:
                vmin_val = np.min(fes_filtered)
                vmax_val = vmax
            else:
                vmin_val = np.min(fes_filtered)
                vmax_val = np.max(fes_filtered)
            
            scatter = ax.scatter(cv1_filtered, cv2_filtered, c=fes_filtered, 
                                cmap=cmap, s=30, alpha=0.8, edgecolors='none',
                                vmin=vmin_val, vmax=vmax_val)
            
            # Set axis ranges
            if cv1_range is not None:
                ax.set_xlim(cv1_range)
            if cv2_range is not None:
                ax.set_ylim(cv2_range)
            
            ax.set_xlabel('CV1', fontsize=14)
            ax.set_ylabel('CV2', fontsize=14)
            ax.set_title(f'Free Energy Surface (γ={self.gamma}, unit: {self.output_unit})', 
                        fontsize=16, pad=20)
            
            cbar = plt.colorbar(scatter, ax=ax)
            cbar.set_label(f'Free Energy ({self.output_unit})', fontsize=14)
            
            # Find minimum (in filtered data)
            if len(fes_filtered) > 0:
                min_idx = np.argmin(fes_filtered)
                min_cv1_val = cv1_filtered[min_idx]
                min_cv2_val = cv2_filtered[min_idx]
                min_fes_val = fes_filtered[min_idx]
            else:
                min_cv1_val = min_cv2_val = min_fes_val = 0
            
            # Generate label based on format
            if min_label_format == 'cv_only':
                min_label = f'Min: ({min_cv1_val:.2f}, {min_cv2_val:.2f})'
            elif min_label_format == 'fes_only':
                min_label = f'Min F: {min_fes_val:.2f} {self.output_unit}'
            elif min_label_format == 'both':
                min_label = f'Min: ({min_cv1_val:.2f}, {min_cv2_val:.2f}), F: {min_fes_val:.2f} {self.output_unit}'
            elif min_label_format == 'cv_fes':
                min_label = f'CV=({min_cv1_val:.2f}, {min_cv2_val:.2f}), F={min_fes_val:.2f} {self.output_unit}'
            else:
                min_label = f'Min: ({min_cv1_val:.2f}, {min_cv2_val:.2f})'
            
            # Mark minimum
            if add_min_marker and len(fes_filtered) > 0:
                ax.plot(min_cv1_val, min_cv2_val, 'r*', markersize=15,
                       markeredgecolor='black', markeredgewidth=1,
                       label=min_label)
                ax.legend(loc='upper right', fontsize=11, framealpha=0.9)
            
            # Mark initial coordinate points
            if initial_coords is not None:
                for i, (init_cv1, init_cv2) in enumerate(initial_coords):
                    # Check if coordinate point is within display range
                    in_x_range = (cv1_range is None) or (cv1_range[0] <= init_cv1 <= cv1_range[1])
                    in_y_range = (cv2_range is None) or (cv2_range[0] <= init_cv2 <= cv2_range[1])
                    
                    if in_x_range and in_y_range:
                        ax.plot(init_cv1, init_cv2, 'go', markersize=10,
                               markeredgecolor='black', markeredgewidth=1,
                               alpha=0.8, label=f'Initial point {i+1}' if i == 0 else "")
            
            # Generate output filename
            if plot_file is None:
                plot_file = 'fes_scatter_heatmap'
            else:
                plot_file = os.path.splitext(plot_file)[0] + '_scatter'
            
            # Save image
            output_formats = output_format.lower().split(',')
            saved_files = []
            
            for fmt in output_formats:
                fmt = fmt.strip()
                output_filename = f"{plot_file}.{fmt}"
                plt.savefig(output_filename, dpi=dpi, bbox_inches='tight')
                saved_files.append(output_filename)
                print(f"Scatter heatmap saved as: {output_filename}")
            
            plt.show()
            
            return fig, ax, saved_files

def parse_range(value):
    """Parse range string, e.g., '0.0:10.0' or '0.0,10.0'"""
    if value is None:
        return None
    
    # Replace colon with comma
    value = value.replace(':', ',')
    
    try:
        # Attempt to parse as list
        values = eval(value)
        if isinstance(values, list) and len(values) == 2:
            return [float(values[0]), float(values[1])]
        else:
            raise ValueError
    except:
        # If failed, try splitting by comma
        parts = value.split(',')
        if len(parts) == 2:
            return [float(parts[0]), float(parts[1])]
        else:
            raise argparse.ArgumentTypeError(f"Invalid range format: {value}, please use 'min,max' format")

def parse_coords(value):
    """Parse coordinate point string, e.g., '1.0,2.0' or '1.0,2.0;3.0,4.0'"""
    if value is None:
        return None
    
    coords = []
    # Split multiple coordinate points by semicolon
    points = value.split(';')
    
    for point in points:
        parts = point.split(',')
        if len(parts) == 2:
            try:
                cv1 = float(parts[0].strip())
                cv2 = float(parts[1].strip())
                coords.append((cv1, cv2))
            except ValueError:
                raise argparse.ArgumentTypeError(f"Invalid coordinate format: {point}")
        else:
            raise argparse.ArgumentTypeError(f"Coordinate point requires two values, got: {point}")
    
    return coords

def main():
    parser = argparse.ArgumentParser(
        description='WT-MTD Free Energy Surface Reconstruction Tool (supports coordinate range restriction, color bar upper limit, and initial point marking)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Coordinate range restriction examples:
  Restrict CV1 display range: --cv1-range "0.0,10.0"
  Restrict CV2 display range: --cv2-range "-5.0,5.0"
  
Color bar upper limit setting:
  Set free energy display upper limit to 20 kcal/mol: --vmax 20.0

Initial coordinate point marking:
  Mark one initial point: --initial-coords "1.5,2.3"
  Mark multiple initial points: --initial-coords "1.5,2.3;3.2,4.1;5.0,1.0"

Unit conversion examples:
  Assume negbias file is in kJ/mol, want kcal/mol output:
    %(prog)s --negbias negbias.dat --gamma 10.0 --output-unit kcal/mol
  
  If negbias file is already in kcal/mol:
    %(prog)s --negbias negbias.dat --gamma 10.0 --input-unit kcal/mol

Contour lines and output format:
  Add contour lines and output EPS:
    %(prog)s --negbias negbias.dat --gamma 10.0 --contour --format eps
  
  Output both PNG and EPS:
    %(prog)s --negbias negbias.dat --gamma 10.0 --format png,eps

Custom contour levels:
    %(prog)s --negbias negbias.dat --gamma 10.0 --contour-levels "[0, 5, 10, 15, 20]"

Minimum label format:
  --min-label cv_only    : Show only CV coordinates (default)
  --min-label fes_only   : Show only free energy value
  --min-label both       : Show CV coordinates and free energy value
  --min-label cv_fes     : Show "CV=(x,y), F=value"

Common unit conversions:
  1 kJ/mol = 0.239006 kcal/mol
  1 kcal/mol = 4.184 kJ/mol
        """
    )
    
    parser.add_argument('--negbias', required=True,
                       help='Path to negbias input file')
    parser.add_argument('--gamma', type=float, required=True,
                       help='Bias factor γ (required)')
    parser.add_argument('--temp', type=float, default=300.0,
                       help='Temperature (K), default 300K')
    parser.add_argument('--kT', type=float, default=None,
                       help='Directly specify kT value (input units)')
    parser.add_argument('--input-unit', choices=['kJ/mol', 'kcal/mol'], default='kJ/mol',
                       help='Unit of input negbias file, default kJ/mol')
    parser.add_argument('--output-unit', choices=['kJ/mol', 'kcal/mol'], default='kcal/mol',
                       help='Unit of output free energy, default kcal/mol')
    parser.add_argument('--output', default=None,
                       help='Output filename (optional)')
    parser.add_argument('--plot', action='store_true',
                       help='Plot heatmap')
    parser.add_argument('--no-shift', action='store_true',
                       help='Do not shift minimum to zero')
    parser.add_argument('--cmap', default='RdYlBu_r',
                       help='Colormap name, default RdYlBu_r')
    parser.add_argument('--dpi', type=int, default=300,
                       help='Output image DPI, default 300')
    parser.add_argument('--no-min-marker', action='store_true',
                       help='Do not mark minimum in heatmap')
    parser.add_argument('--grid', action='store_true',
                       help='Show grid in heatmap')
    parser.add_argument('--min-label', choices=['cv_only', 'fes_only', 'both', 'cv_fes'], 
                       default='cv_only',
                       help='Minimum label format, default cv_only')
    parser.add_argument('--contour', action='store_true',
                       help='Add contour lines')
    parser.add_argument('--contour-levels', type=str, default='10',
                       help='Number of contour levels or list, e.g., "10" or "[0,5,10,15,20]", default 10')
    parser.add_argument('--contour-colors', default='black',
                       help='Contour line color, default black')
    parser.add_argument('--contour-alpha', type=float, default=0.7,
                       help='Contour line transparency, default 0.7')
    parser.add_argument('--contour-linewidths', type=float, default=0.8,
                       help='Contour line width, default 0.8')
    parser.add_argument('--format', default='png',
                       help='Output image format, multiple formats separated by commas, e.g., "png,eps,pdf", default png')
    
    # New parameters
    parser.add_argument('--cv1-range', type=parse_range, default=None,
                       help='CV1 display range, format: min,max or min:max, e.g., "0.0,10.0"')
    parser.add_argument('--cv2-range', type=parse_range, default=None,
                       help='CV2 display range, format: min,max or min:max, e.g., "-5.0,5.0"')
    parser.add_argument('--vmax', type=float, default=None,
                       help='Color bar maximum value (free energy upper limit), e.g., 20.0 displays 0-20 kcal/mol range')
    parser.add_argument('--initial-coords', type=parse_coords, default=None,
                       help='Initial coordinate points, format: "cv1,cv2" or "cv1,cv2;cv1,cv2;...", e.g., "1.5,2.3;3.2,4.1"')
    
    args = parser.parse_args()
    
    # Check if file exists
    if not os.path.exists(args.negbias):
        print(f"Error: File does not exist: {args.negbias}")
        return
    
    # Parse contour levels
    try:
        if args.contour_levels.startswith('[') and args.contour_levels.endswith(']'):
            # Parse list
            contour_levels = eval(args.contour_levels)
        else:
            # Parse integer
            contour_levels = int(args.contour_levels)
    except:
        print(f"Warning: Unable to parse contour levels '{args.contour_levels}', using default value 10")
        contour_levels = 10
    
    # Create reconstructor
    reconstructor = WT_FES_Reconstructor(
        gamma=args.gamma,
        input_unit=args.input_unit,
        output_unit=args.output_unit,
        kT=args.kT,
        temp=args.temp
    )
    
    # Reconstruct free energy surface
    cv1, cv2, fes, output_file = reconstructor.reconstruct_from_negbias(
        negbias_file=args.negbias,
        output_file=args.output,
        shift_min_to_zero=not args.no_shift
    )
    
    # Plotting
    saved_files = []
    if args.plot:
        plot_file = os.path.splitext(output_file)[0] + '.png'
        fig, ax, saved_files = reconstructor.plot_fes_heatmap(
            cv1, cv2, fes, 
            plot_file=plot_file,
            cmap=args.cmap,
            dpi=args.dpi,
            add_min_marker=not args.no_min_marker,
            add_grid=args.grid,
            min_label_format=args.min_label,
            add_contour=args.contour,
            contour_levels=contour_levels,
            contour_colors=args.contour_colors,
            contour_alpha=args.contour_alpha,
            contour_linewidths=args.contour_linewidths,
            output_format=args.format,
            # New parameters
            cv1_range=args.cv1_range,
            cv2_range=args.cv2_range,
            vmax=args.vmax,
            initial_coords=args.initial_coords
        )
    
    print(f"\n{'='*60}")
    print("Reconstruction completed!")
    print(f"Output file: {output_file}")
    if args.plot:
        print(f"Image files:")
        for file in saved_files:
            print(f"  {file}")
    
    # Print key information
    print(f"\nKey conversion information:")
    print(f"  Input file unit: {args.input_unit}")
    print(f"  Output free energy unit: {args.output_unit}")
    print(f"  Bias factor γ = {args.gamma}")
    print(f"  WT-MTD conversion factor = {reconstructor.conversion_factor:.6f}")
    print(f"  Unit conversion factor = {reconstructor.unit_conversion_factor:.6f}")
    
    # New feature information
    if args.cv1_range is not None:
        print(f"  CV1 display range: [{args.cv1_range[0]}, {args.cv1_range[1]}]")
    if args.cv2_range is not None:
        print(f"  CV2 display range: [{args.cv2_range[0]}, {args.cv2_range[1]}]")
    if args.vmax is not None:
        print(f"  Color bar maximum: {args.vmax} {args.output_unit}")
    if args.initial_coords is not None:
        print(f"  Number of initial coordinate points: {len(args.initial_coords)}")
    
    print(f"  Contour lines: {'Enabled' if args.contour else 'Disabled'}")
    if args.contour:
        print(f"  Contour levels: {contour_levels}")
    print(f"  Output format: {args.format}")
    print(f"  Minimum label format: {args.min_label}")
    
    # Unit conversion hints
    if args.input_unit != args.output_unit:
        print(f"\nUnit conversion:")
        if args.input_unit == 'kJ/mol' and args.output_unit == 'kcal/mol':
            print(f"  1 kJ/mol = 0.239006 kcal/mol")
            print(f"  Typical protein folding barrier: ~50 kJ/mol ≈ 11.95 kcal/mol")
        else:
            print(f"  1 kcal/mol = 4.184 kJ/mol")
            print(f"  Typical protein folding barrier: ~12 kcal/mol ≈ 50.21 kJ/mol")

if __name__ == "__main__":
    main()