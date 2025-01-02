import json

# Function to detect peaks in a spectrum
def detect_peaks(spectrum):
    peaks = []
    for i in range(1, len(spectrum) - 1):
        if float(spectrum[i]["Rel. Intensity"]) > float(spectrum[i - 1]["Rel. Intensity"]) and \
                float(spectrum[i]["Rel. Intensity"]) > float(spectrum[i + 1]["Rel. Intensity"]):
            peaks.append(i)
    return peaks


# Function to calculate peak areas (intensity sums) in a spectrum
def calculate_peak_areas(spectrum, peaks):
    peak_areas = []
    for i in range(len(peaks) - 1):
        start_index = peaks[i]
        end_index = peaks[i + 1]
        area = sum(float(entry["Rel. Intensity"]) for entry in spectrum[start_index:end_index])
        peak_areas.append(area)
    return peak_areas

# Load the JSON data
with open('tokenised-full-27-2.json', 'r') as file:
    mass_spec_data = json.load(file)


# Consolidated mass spec data list
consolidated_mass_spec_data = []

# Iterate through the mass spec data
current_entry = None
for entry in mass_spec_data:
    if current_entry is None:
        current_entry = entry
        continue
    
    # Check if the current entry and the previous entry have the same OriginalMolecule value
    if entry['OriginalMolecule'] == current_entry['OriginalMolecule']:
        # Merge spectra
        current_entry['Spectrum'] += entry['Spectrum']
    else:
        # Normalize intensity values based on peak areas (sums)
        peaks = detect_peaks(current_entry['Spectrum'])
        peak_areas = calculate_peak_areas(current_entry['Spectrum'], peaks)
        total_peak_area = sum(peak_areas)
        
        print("Type of current_entry['Spectrum']:", type(current_entry['Spectrum']))
        print("Type of total_peak_area:", type(total_peak_area))

        normalized_spectrum = [value / total_peak_area for value in current_entry['Spectrum']]
        current_entry['Spectrum'] = normalized_spectrum
        
        # Add the previous entry to the consolidated list
        consolidated_mass_spec_data.append(current_entry)
        current_entry = entry

# Add the last entry to the consolidated list
if current_entry is not None:
    # Normalize intensity values based on peak areas (sums)
    peaks = detect_peaks(current_entry['Spectrum'])
    peak_areas = calculate_peak_areas(current_entry['Spectrum'], peaks)
    total_peak_area = sum(peak_areas)
    normalized_spectrum = [value / total_peak_area for value in current_entry['Spectrum']]
    current_entry['Spectrum'] = normalized_spectrum
    
    consolidated_mass_spec_data.append(current_entry)

# Calculate the total number of entries in the consolidated list
total_entries = len(consolidated_mass_spec_data)

# Write the consolidated mass spec data to a new JSON file
with open('consolidated_normalized_mass_spec_data.json', 'w') as file:
    json.dump(consolidated_mass_spec_data, file, indent=4)

print("Total number of entries in the consolidated list:", total_entries)
