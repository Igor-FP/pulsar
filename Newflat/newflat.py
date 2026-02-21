#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
newflat.py - Add maintenance/flat renewal entry to technical log.

Records the timestamp when camera optics were changed, requiring new flats.
Used by autocalibrate.py to determine flat validity periods.

Usage:
    newflat.py --camera "2600MM" --log maintenance.csv
    newflat.py --camera "2600MM" --log maintenance.csv --date "2024-05-18T14:30:00"

Log format (CSV):
    2024-05-18T14:30:00,2600MM
    2024-12-21T10:00:00,2600MM

The camera string is matched as substring against INSTRUME header field.
"""

import sys
import os
import argparse
from datetime import datetime, timezone


def parse_args():
    parser = argparse.ArgumentParser(
        description="Add maintenance entry to flat validity log",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    newflat.py --camera "2600MM" --log maintenance.csv
        Add entry with current UTC time

    newflat.py --camera "2600MM" --log maintenance.csv --date "2024-05-18T14:30:00"
        Add entry with specific UTC time

    newflat.py --camera "QHY600" --log /path/to/site_a.csv
        Different camera/log for another setup

Log format:
    Each line: DATETIME_UTC,CAMERA_ID
    Example: 2024-05-18T14:30:00,2600MM
"""
    )

    parser.add_argument(
        "--camera", "-c",
        required=True,
        help="Camera identifier (substring to match INSTRUME header, e.g. '2600MM')"
    )

    parser.add_argument(
        "--log", "-l",
        required=True,
        help="Path to maintenance log CSV file (created if not exists)"
    )

    parser.add_argument(
        "--date", "-d",
        default=None,
        help="UTC datetime (ISO format: YYYY-MM-DDTHH:MM:SS). Default: current UTC time"
    )

    parser.add_argument(
        "--comment", "-m",
        default=None,
        help="Optional comment (will be added as third column)"
    )

    return parser.parse_args()


def main():
    args = parse_args()

    # Determine timestamp
    if args.date:
        try:
            dt = datetime.fromisoformat(args.date.replace("Z", "+00:00"))
            # If no timezone specified, assume UTC
            if dt.tzinfo is None:
                dt = dt.replace(tzinfo=timezone.utc)
        except ValueError as e:
            print(f"ERROR: Invalid date format: {args.date}")
            print("Expected ISO format: YYYY-MM-DDTHH:MM:SS")
            sys.exit(1)
    else:
        dt = datetime.now(timezone.utc)

    # Format without timezone suffix (UTC implied)
    timestamp = dt.strftime("%Y-%m-%dT%H:%M:%S")

    # Sanitize camera ID (no commas or newlines)
    camera = args.camera.replace(",", "_").replace("\n", " ").strip()
    if not camera:
        print("ERROR: Camera identifier cannot be empty")
        sys.exit(1)

    # Build log line
    if args.comment:
        comment = args.comment.replace(",", ";").replace("\n", " ").strip()
        line = f"{timestamp},{camera},{comment}\n"
    else:
        line = f"{timestamp},{camera}\n"

    # Ensure directory exists
    log_dir = os.path.dirname(args.log)
    if log_dir and not os.path.exists(log_dir):
        os.makedirs(log_dir, exist_ok=True)

    # Check if file exists (for header)
    file_exists = os.path.isfile(args.log)

    # Append to log
    try:
        with open(args.log, "a", encoding="utf-8") as f:
            # Add header comment if new file
            if not file_exists:
                f.write("# Flat maintenance log - UTC timestamps\n")
                f.write("# Format: DATETIME_UTC,CAMERA_ID[,COMMENT]\n")
            f.write(line)
    except IOError as e:
        print(f"ERROR: Cannot write to log file: {e}")
        sys.exit(1)

    print(f"Added: {timestamp} | Camera: {camera}")
    if args.comment:
        print(f"Comment: {args.comment}")
    print(f"Log: {os.path.abspath(args.log)}")


if __name__ == "__main__":
    main()
