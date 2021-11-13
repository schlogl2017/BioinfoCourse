#!/usr/bin/env bash

IFS= read -r header        
printf '%s\n' "$header"    
$@   
