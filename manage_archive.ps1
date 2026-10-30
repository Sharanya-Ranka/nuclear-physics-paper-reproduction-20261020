# --- CONFIGURATION ---
# Define the name of the archive file
$archiveName = "data_archive.tar.gz"

# Define the specific files you want to target (paths can be absolute or relative)
$filesToArchive = @(
    "./data_sourcing/plots",
    "./data_sourcing/data_processed",
    "./data_sourcing/data",
    "./pinn/data"
)

# Define the destination directory for extraction
$extractDestination = "."
# ---------------------

function Protect-Archive {
    Write-Host "Creating archive: $archiveName..." -ForegroundColor Cyan
    
    # Use native Windows tar to compress the specific files
    # -c: Create, -z: Compress (gzip), -f: File name
    tar -czf $archiveName $filesToArchive
    
    if ($LASTEXITCODE -eq 0) {
        Write-Host "Successfully created $archiveName!" -ForegroundColor Green
    } else {
        Write-Error "Failed to create archive."
    }
}

function Expand-Archive {
    Write-Host "Extracting archive: $archiveName to $extractDestination..." -ForegroundColor Cyan
    
    # Ensure destination directory exists
    if (!(Test-Path -Path $extractDestination)) {
        New-Item -ItemType Directory -Path $extractDestination | Out-Null
    }
    
    # Use native Windows tar to extract specific files
    # -x: Extract, -z: Uncompress (gzip), -f: File name, -C: Target directory
    tar -xzf $archiveName -C $extractDestination
    
    if ($LASTEXITCODE -eq 0) {
        Write-Host "Successfully extracted files to $extractDestination!" -ForegroundColor Green
    } else {
        Write-Error "Failed to extract archive."
    }
}

# --- EXECUTION ---
# Uncomment the action you want to perform:

Protect-Archive
# Expand-Archive
