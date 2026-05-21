#!/usr/bin/env bash

set -u

# Copy NIfTI files into matching raw fMRI sub/ses folders.
# Default source and target can be overridden with CLI flags.
NIFTI_SOURCE="/Volumes/TOSHIBA_EXT/03-Raw_fMRI/2026/nifti"
RAW_FMRI_TARGET="/Volumes/TOSHIBA_EXT/03-Raw_fMRI"
DRY_RUN=0
FORCE_OVERWRITE=0

GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

print_usage() {
    cat <<EOF
Usage: $0 [options]

Options:
  --source PATH      Source nifti root (default: $NIFTI_SOURCE)
  --target PATH      Raw fMRI target root (default: $RAW_FMRI_TARGET)
  --dry-run          Build and print plan only; do not copy
  --overwrite        Allow overwriting existing destination files
  -h, --help         Show this help message
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --source)
            NIFTI_SOURCE="$2"
            shift 2
            ;;
        --target)
            RAW_FMRI_TARGET="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        --overwrite)
            FORCE_OVERWRITE=1
            shift
            ;;
        -h|--help)
            print_usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            print_usage
            exit 2
            ;;
    esac
done

if [[ ! -d "$NIFTI_SOURCE" ]]; then
    echo -e "${RED}ERROR:${NC} Source path not found: $NIFTI_SOURCE"
    exit 1
fi

if [[ ! -d "$RAW_FMRI_TARGET" ]]; then
    echo -e "${RED}ERROR:${NC} Target path not found: $RAW_FMRI_TARGET"
    exit 1
fi

categorize_file() {
    local filename="$1"
    local lower
    lower=$(printf '%s' "$filename" | tr '[:upper:]' '[:lower:]')

    # Fieldmap first so e1/e2/ph names are not misclassified as func.
    if [[ "$lower" == *"fieldmap"* ]] || [[ "$lower" == *"field_mapping"* ]] || [[ "$lower" == *"gre field mapping"* ]] || [[ "$lower" == *"gre_field_mapping"* ]] || [[ "$lower" == *"fmri_fieldmap"* ]] || [[ "$lower" == *"_e1"* ]] || [[ "$lower" == *"_e2"* ]] || [[ "$lower" == *"_ph"* ]] || [[ "$lower" == *"phase"* ]] || [[ "$lower" == *"magnitude"* ]]; then
        printf 'fmap\n'
        return
    fi

    if [[ "$lower" == *"3d_t1"* ]] || [[ "$lower" == *"3d t1"* ]] || [[ "$lower" == *"t1w"* ]] || [[ "$lower" == *"mprage"* ]] || [[ "$lower" == *"anat"* ]]; then
        printf 'anat\n'
        return
    fi

    if [[ "$lower" == *"fmri"* ]] || [[ "$lower" == *"memorytask"* ]] || [[ "$lower" == *"memory_task"* ]] || [[ "$lower" == *"task"* ]] || [[ "$lower" == *"bold"* ]] || [[ "$lower" == *"rest"* ]]; then
        printf 'func\n'
        return
    fi

    printf 'unknown\n'
}

copy_file_with_permission_fallback() {
    local src="$1"
    local dst_dir="$2"
    local err_file="$3"
    local cp_args=(-p)

    if [[ "$FORCE_OVERWRITE" -eq 0 ]]; then
        cp_args+=(-n)
    fi

    if cp "${cp_args[@]}" "$src" "$dst_dir/" 2>"$err_file"; then
        return 0
    fi

    if grep -Eiq 'permission denied|operation not permitted|not permitted' "$err_file"; then
        if sudo cp "${cp_args[@]}" "$src" "$dst_dir/" 2>"$err_file"; then
            return 0
        fi
    fi

    return 1
}

COPY_PLAN_FILE=$(mktemp)
FAILED_FILE=$(mktemp)
trap 'rm -f "$COPY_PLAN_FILE" "$FAILED_FILE"' EXIT

echo "============================================"
echo "NIfTI to Raw fMRI Copy Script"
echo "============================================"
echo "Source: $NIFTI_SOURCE"
echo "Target: $RAW_FMRI_TARGET"
echo "Dry run: $DRY_RUN"
echo "Overwrite: $FORCE_OVERWRITE"
echo

echo -e "${YELLOW}Step 1/3:${NC} Building copy plan"

planned=0
skipped_existing=0
skipped_no_target=0
skipped_unknown=0
source_folders_seen=0

while IFS= read -r -d '' nifti_folder; do
    source_folders_seen=$((source_folders_seen + 1))
    folder_name=$(basename "$nifti_folder")

    if [[ "$folder_name" =~ ^([0-9]+)_irm_(t[0-9]+)$ ]]; then
        subject_id="${BASH_REMATCH[1]}"
        session="${BASH_REMATCH[2]}"
    else
        echo -e "${YELLOW}WARN:${NC} Skipping unmatched source folder name: $folder_name"
        continue
    fi

    target_base="$RAW_FMRI_TARGET/sub-$subject_id/ses-$session"
    if [[ ! -d "$target_base" ]]; then
        skipped_no_target=$((skipped_no_target + 1))
        echo -e "${YELLOW}WARN:${NC} Missing target folder for $folder_name -> $target_base"
        continue
    fi

    while IFS= read -r -d '' source_file; do
        filename=$(basename "$source_file")
        category=$(categorize_file "$filename")

        if [[ "$category" == "unknown" ]]; then
            skipped_unknown=$((skipped_unknown + 1))
            echo -e "${YELLOW}WARN:${NC} Unknown category, skipping file: $source_file"
            continue
        fi

        target_dir="$target_base/$category"
        target_file="$target_dir/$filename"

        if [[ "$FORCE_OVERWRITE" -eq 0 && -e "$target_file" ]]; then
            skipped_existing=$((skipped_existing + 1))
            continue
        fi

        printf '%s\t%s\t%s\n' "$source_file" "$target_dir" "$category" >> "$COPY_PLAN_FILE"
        planned=$((planned + 1))
    done < <(find "$nifti_folder" -maxdepth 1 -type f \( -iname '*.nii' -o -iname '*.nii.gz' \) -print0)
done < <(find "$NIFTI_SOURCE" -mindepth 1 -maxdepth 1 -type d -name '*_irm_t*' -print0)

echo -e "${GREEN}Plan complete.${NC}"
echo "  Source folders scanned: $source_folders_seen"
echo "  Files planned: $planned"
echo "  Skipped (existing in target): $skipped_existing"
echo "  Skipped (no matching target folder): $skipped_no_target"
echo "  Skipped (unknown category): $skipped_unknown"
echo

if [[ "$planned" -eq 0 ]]; then
    echo -e "${YELLOW}Nothing to copy.${NC}"
    exit 0
fi

echo -e "${BLUE}Planned operations:${NC}"
awk -F '\t' '{print "  - " $1 " -> " $2 "/"}' "$COPY_PLAN_FILE"
echo

if [[ "$DRY_RUN" -eq 1 ]]; then
    echo -e "${GREEN}Dry run complete.${NC}"
    exit 0
fi

echo -e "${YELLOW}Step 2/3:${NC} Executing copy plan"

copied=0
failed=0
index=0

while IFS=$'\t' read -r source_file target_dir category; do
    index=$((index + 1))
    filename=$(basename "$source_file")
    mkdir_err=$(mktemp)

    if ! mkdir -p "$target_dir" 2>"$mkdir_err"; then
        if grep -Eiq 'permission denied|operation not permitted|not permitted' "$mkdir_err"; then
            if ! sudo mkdir -p "$target_dir" 2>"$mkdir_err"; then
                failed=$((failed + 1))
                printf '%s\t%s\t%s\n' "$source_file" "$target_dir" "$(cat "$mkdir_err")" >> "$FAILED_FILE"
                rm -f "$mkdir_err"
                echo -e "[$index/$planned] ${RED}FAIL${NC} mkdir: $target_dir"
                continue
            fi
        else
            failed=$((failed + 1))
            printf '%s\t%s\t%s\n' "$source_file" "$target_dir" "$(cat "$mkdir_err")" >> "$FAILED_FILE"
            rm -f "$mkdir_err"
            echo -e "[$index/$planned] ${RED}FAIL${NC} mkdir: $target_dir"
            continue
        fi
    fi
    rm -f "$mkdir_err"

    cp_err=$(mktemp)
    if copy_file_with_permission_fallback "$source_file" "$target_dir" "$cp_err"; then
        copied=$((copied + 1))
        echo -e "[$index/$planned] ${GREEN}OK${NC} $filename -> $target_dir/ ($category)"
    else
        failed=$((failed + 1))
        printf '%s\t%s\t%s\n' "$source_file" "$target_dir" "$(cat "$cp_err")" >> "$FAILED_FILE"
        echo -e "[$index/$planned] ${RED}FAIL${NC} $filename -> $target_dir/"
    fi
    rm -f "$cp_err"
done < "$COPY_PLAN_FILE"

echo
echo -e "${YELLOW}Step 3/3:${NC} Final summary"
echo "============================================"
echo "Planned files: $planned"
echo -e "Copied: ${GREEN}$copied${NC}"
echo -e "Failed: ${RED}$failed${NC}"

if [[ "$failed" -gt 0 ]]; then
    echo
    echo "Failures (source -> target -> error):"
    awk -F '\t' '{print "  - " $1 " -> " $2 " -> " $3}' "$FAILED_FILE"
    exit 1
fi

echo -e "${GREEN}All planned files copied successfully.${NC}"
exit 0
