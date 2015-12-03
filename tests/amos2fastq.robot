*** Settings ***
Library         Process
Library         OperatingSystem
Library         Collections 
Suite Teardown          Terminate All Processes    

*** Variables ***
@{T_EXPECTED}       expected_fooCTG0.fastq           expected_fooCTG1.fastq         expected_fooCTG2.fastq          expected_fooCTG3.fastq          expected_fooCTG4.fastq
@{T_ACTUAL}       foo-0.fastq          foo-1.fastq        foo-2.fastq           foo-3.fastq      foo-4.fastq
${in_fastq1} =          tests/testinput/foo1.fastq
${in_fastq2} =          tests/testinput/foo2.fastq
${in_amos} =            tests/testinput/foo.afg

*** Test Cases ***
TestAmos2Fastq
    ${process_result} =         Run Process     python          bio_bits/amos2fastq_main.py      ${in_fastq1}            ${in_fastq2}    --amos      ${in_amos}
    # Check system exited  correctly
    Should Be Equal As Integers         ${process_result.rc}        0 
    Log To Console       ${process_result.stdout}
    Log To Console        ${process_result.stderr}
    # Check output
    Should Not Contain         ${process_result.stdout}        Error

# This can become generic enough to load variables from a file (expected and actual filenames)
# Then assert equality for all of them.
# pybot is NOT case sensitive.

        : FOR       ${i}          IN RANGE        4
        \       Log to Console          ${i}
        \        ${expected} =       Get From List   ${T_EXPECTED}       ${i}
        \        ${actual} =     Get From List   ${T_ACTUAL}       ${i}
        \        Log To Console    ${expected}
        \        Log To Console    ${actual}
        \        File Should Exist       ${actual}
        \        File Should Not Be Empty        ${actual}
        \        ${actual_contents} =    Get File        ${actual}
        \        ${expected_contents} =  Get File        tests/expected/${expected}
        \        Should Be Equal As Strings     ${actual_contents}      ${expected_contents}
        \       Log To Console          passed
        
