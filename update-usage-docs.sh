cp doc/source/usage_template.rst usage_template.rst

for program in projectmanager trackmanager; do
    for cmd in "--help"; do
        echo $program
        echo $cmd
        echo "\n\`\`$program $cmd\`\`" > ${program}_USAGE_header.temp
        # echo -e "----------------------------------" >> ${program}_USAGE_header.temp
        $program $cmd --help > ${program}_USAGE.temp 2>&1
        sed -i 's/^/\t/' ${program}_USAGE.temp
        sed -i '1s/^/\n.. code-block:: none\n\n/' ${program}_USAGE.temp
        #sed -i -e "/\`looper $cmd\`/r ${program}_USAGE.temp" -e '$G' ${program}_usage_template.rst  # for -in place inserts
        cat ${program}_USAGE_header.temp ${program}_USAGE.temp >> usage_template.rst # to append to the end
    done
    rm ${program}_USAGE.temp
    rm ${program}_USAGE_header.temp
done

mv usage_template.rst  doc/source/manager_programs.rst
