*** Settings ***
Library         Process
Library         OperatingSystem
Library         Collections 
Suite Teardown          Terminate All Processes    

*** Variables ***
${EXPECTED} =      tests/testinput/ctl_expected.tsv
${ACTUAL} =     robotout.tsv
${in_fasta} =   tests/Den4_MAAPS_TestData16.fasta 
${in_genbank} =         tests/testinput/sequence.gb

*** Test Cases ***
Testdegen_regions
    ${process_result} =         Run Process     degen_regions        -i      ${in_fasta}     -o      ${ACTUAL}       --gb-file       ${in_genbank}

    # Check system exited  correctly
    Should Be Equal As Integers         ${process_result.rc}        0 
    Log To Console       ${process_result.stdout}
    Log To Console        ${process_result.stderr}

    # Check output
    Should Not Contain         ${process_result.stdout}        Error
    File Should Exist       ${actual}
    File Should Not Be Empty        ${actual}
    ${actual_contents} =    Get File        ${actual}
    ${expected_contents} =  Get File        ${expected}
    Should Be Equal As Strings     ${actual_contents}      ${expected_contents} 
